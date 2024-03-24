
# const TRACE_YRANGES = Dict{String, Vector{Int64}}(
#     "J" => [756, 1690],
#     "Y" => [1020, 1780],
#     "H" => [85, 2020],
#     #"JH" => [], # manual
#     # YJ => [], # manual
#     "KCEN" => [58, 558],
#     "HK" => [121, 1821],
# )


const DEFAULT_REDUCE_OPTIONS = Dict(
    "mode" => nothing,
    "trace_xrange" => [10, 2048-10],
    "trace_yrange" => nothing,
    "trace_offset" => 0,
    "trace_fit_window" => 100,
    "trace_poly_deg" => 2,
    "trace_width_bounds" => [30, 120],
    "trace_peak_thresh" => 0.3,
    "trace_σ_bounds" => [1, 10],
    "trace_x0" => 1024,
    "trace_flux_smooth_width" => 800,
    "trace_min_spacing" => 2,
    "extract_xrange" => [200, 2048-200],
    "extract_orders" => nothing,
    "extract_badpix_nσ" => 10,
    "do_flat" => true
)


function reduce(;
        science::Vector{NIRSPECL0},
        darks::Union{Vector{NIRSPECL0}, Nothing}=nothing,
        flats::Union{Vector{NIRSPECL0}, Nothing}=nothing,
        traces::Union{OrderedDict, Nothing}=nothing,
        output_path::String,
        options::Union{AbstractDict{String, <:Any}, Nothing}=nothing
    )

    # Get options
    options = isnothing(options) ? copy(DEFAULT_REDUCE_OPTIONS) : merge(DEFAULT_REDUCE_OPTIONS, options)

    # Dirs
    if output_path[end] != PATHSEP
        output_path *= PATHSEP
    end
    sci_output_path = output_path * "spectra" * PATHSEP
    cal_output_path = output_path * "calib" * PATHSEP
    mkpath(sci_output_path)
    mkpath(cal_output_path)

    # Gen calib images
    median_dark = !isnothing(darks) ? gen_median_dark(darks, cal_output_path, options) : nothing
    median_flat = !isnothing(flats) ? gen_median_flat(flats, cal_output_path, options, median_dark) : nothing

    # Trace the orders
    if isnothing(traces)
        traces = trace(median_flat, cal_output_path, options)
    end

    # Fix regions between the orders
    # Remove blaze from flats
    if options["do_flat"]
        med_flat, continua = remove_slitflat_continuum(median_flat, traces, cal_output_path, options)
    else
        med_flat = nothing
        continua = nothing
    end

    # Extract
    extract(science, traces, sci_output_path, options; dark=median_dark, flat=med_flat, continua=continua)

end


function EchelleReduce.gen_median_dark(darks::Vector{NIRSPECL0}, output_path::String, options::Dict{String, <:Any})
    median_dark = NIRSPECL0(get_median_dark_filename(darks, output_path))
    println("Generating median dark: $median_dark")
    dark_images = read.(darks)
    median_dark_image = EchelleReduce.gen_median_dark(dark_images)
    FITSIO.fitswrite(median_dark.filename, median_dark_image, header=deepcopy(read_header(darks[1])))
    return median_dark
end


function EchelleReduce.gen_median_flat(flats::Vector{NIRSPECL0}, output_path::String, options::Dict{String, <:Any}, dark::Union{NIRSPECL0, Nothing}=nothing)
    median_flat = NIRSPECL0(get_median_flat_filename(flats, output_path))
    println("Generating median flat: $median_flat")
    flat_images = read.(flats)
    dark_image = !isnothing(dark) ? read(dark) : nothing
    flat_itime = get_itime(flats[1])
    dark_itime = !isnothing(dark) ? get_itime(dark) : nothing
    median_flat_image = gen_median_flat(flat_images; dark_image, q=0.75, flat_itime, dark_itime)
    FITSIO.fitswrite(median_flat.filename, median_flat_image, header=deepcopy(read_header(flats[1])))
    return median_flat
end

function get_median_flat_filename(flats::Vector{NIRSPECL0}, output_path::String)
    img_nums = get_image_num.(flats)
    img_start, img_end = minimum(img_nums), maximum(img_nums)
    filename = "$(output_path)median_flat_$(get_utdate(flats[1]))_imgs$(img_start)-$(img_end).fits"
    return filename
end


function get_median_dark_filename(darks::Vector{NIRSPECL0}, output_path::String)
    img_nums = get_image_num.(darks)
    img_start, img_end = minimum(img_nums), maximum(img_nums)
    itime = get_itime(darks[1])
    filename = "$(output_path)median_dark_$(get_utdate(darks[1]))_imgs$(img_start)-$(img_end)_$(itime)s.fits"
    return filename
end


function trace(
        flat::NIRSPECL0,
        output_path::String, options::AbstractDict{String, <:Any},
    )

    println("Tracing Orders: $flat")

    # Load image
    flat_image = read(flat)

    # Labels
    labels = string.(options["trace_orders"])

    # Trace
    traces = trace_boxcar(flat_image, labels;
                    xrange=options["trace_xrange"], yrange=options["trace_yrange"], x0=options["trace_x0"],
                    width_bounds=options["trace_width_bounds"], σ_bounds=options["trace_σ_bounds"],
                    min_spacing=options["trace_min_spacing"], deg=options["trace_poly_deg"], fit_window=options["trace_fit_window"], flux_smooth_width=options["trace_flux_smooth_width"],
                    peak_thresh=options["trace_peak_thresh"]
                )

    # Plot image
    filename_out = "$(output_path)$(basename(flat)[1:end-5])_tracepos.png"
    plot_tracepos_image(flat_image, getproperty.(values(traces), :yc), filename_out, qmax=0.9)

    # Save trace params to jld file
    filename_out = "$(output_path)$(splitext(basename(flat))[1])_traces.jld"
    jldsave(filename_out; traces)

    # Return
    return traces
end


function extract(
        data::NIRSPECL0, traces::OrderedDict{String, <:Any}, options::Dict{String, <:Any};
        dark::Union{NIRSPECL0, Nothing}=nothing,
        flat::Union{NIRSPECL0, Nothing}=nothing,
        continua::Union{<:AbstractDict{String, <:Any}, Nothing}=nothing
    )

    # Load median dark and flat
    dark_image = !isnothing(dark) ? read(dark) : nothing
    flat_image = !isnothing(flat) && options["do_flat"] ? read(flat) : nothing

    # Load data
    data_image = read(data)

    # Calibrate the image
    itime = get_itime(data)
    dark_itime = !isnothing(dark) ? get_itime(dark) : nothing
    calibrate_image!(data_image; dark_image, itime, dark_itime, flat_image)

    # Fix negatives
    data_image[data_image .< 0] .= NaN

    # Convert to pe
    data_image .*= DETECTOR_GAIN

    # Extract traces
    reduced_data = OrderedDict{String, Any}()
    for trace in values(traces)
        if parse(Int, trace.label) in options["extract_orders"]
            extractor = get_extractor(data, traces[trace.label], options)
            reduced_data[trace.label] = extract_trace(data.filename, data_image, trace.label, extractor)
        end
    end

    # Correct continuum
    # if options["do_flat"] && !isnothing(continua)
    #     for trace in values(traces)
    #         if parse(Int, trace.label) in options["extract_orders"] && !isnothing(reduced_data[trace.label])
    #             reduced_data[trace.label].spec ./= continua[trace.label]
    #             reduced_data[trace.label].specerr ./= continua[trace.label]
    #         end
    #     end
    # end

    # Return
    return reduced_data

end


function extract(
        science::Vector{NIRSPECL0},
        traces::OrderedDict{String, <:Any},
        output_path::String, options::Dict{String, <:Any};
        dark::Union{NIRSPECL0, Nothing}=nothing,
        flat::Union{NIRSPECL0, Nothing}=nothing,
        continua::Union{AbstractDict{String, <:Any}, Nothing}=nothing,
    )
    pmap(science) do sci_data
        reduced_data = extract(sci_data, traces, options; dark, flat, continua)
        save_extraction_results(sci_data, reduced_data, output_path)
    end
end


function save_extraction_results(data::NIRSPECL0, reduced_data::OrderedDict{String, <:Any}, output_path::String)

    # Plot
    target = replace(get_object(data), " " => "_")
    filename_out = "$(output_path)$(splitext(basename(data))[1])_$(target)_reduced.png"
    plot_reduced_spectrum(reduced_data, filename_out)

    # Collect extracted spectra
    nx = 2048
    reduced_data_out = OrderedDict{String, Any}()
    for trace_label ∈ keys(reduced_data)
        if !isnothing(reduced_data[trace_label])
            reduced_data_out[trace_label] = [reduced_data[trace_label].spec reduced_data[trace_label].specerr]
        else
            reduced_data_out[trace_label] = fill(NaN, (nx, 2))
        end
    end

    # Save .fits file
    target = replace(get_object(data), " " => "_")
    fname = "$(output_path)$(splitext(basename(data))[1])_$(target)_reduced.fits"
    header = read_header(data)
    FITSIO.FITS(fname, "w") do f
        write(f, reduced_data_out, header=header)
    end

    # Save auxiliary data
    fname = "$(output_path)$(splitext(basename(data))[1])_$(target)_auxiliary.jld"
    save_auxiliary_data(fname, reduced_data)

end


function get_read_noise(::AnyNIRSPEC)
    return DETECTOR_READ_NOISE
end


function EchelleReduce.remove_slitflat_continuum(flat::NIRSPECL0, traces::OrderedDict{String, <:NamedTuple}, output_path::String, options::Dict{String, <:Any})

    # Load header/image
    header = read_header(flat)
    flat_image = read(flat)

    # Remove continuum
    flat_image_norm, continua = remove_slitflat_continuum(flat_image, traces; xrange=options["trace_xrange"], deg=4)

    # Save new flat
    flat_norm = NIRSPECL0(output_path * basename(flat.filename)[1:end-5] * "_norm.fits")
    FITSIO.fitswrite(
        flat_norm.filename,
        flat_image_norm,
        header=header
    )
    jldsave(output_path * basename(flat.filename)[1:end-5] * "_continua.jld"; continua)

    # Return
    return flat_norm, continua

end


function get_extractor(data::NIRSPECL0, trace::NamedTuple, options::Dict{String, <:Any})
    read_noise = get_read_noise(data)
    return Boxcar(;
        xrange=options["extract_xrange"],
        yrange=options["extract_yrange"],
        yc=trace.yc .+ options["trace_offset"],
        collapse_function="median",
        read_noise=read_noise,
    )
end

