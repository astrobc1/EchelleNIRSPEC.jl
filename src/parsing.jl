

# Read headers
Echelle.read_header(data::NIRSPECL0) = read_header(data.filename, 1)
Echelle.read_header(data::NIRSPECL1) = read_header(data.filename, 2)
Echelle.read_key(data::NIRSPECL0, key::Union{Int, String}) = read_key(data.filename, key, hdu=1)
Echelle.read_key(data::NIRSPECL1, key::Union{Int, String}) = read_key(data.filename, key, hdu=2)


# Parsing
Echelle.get_object(data::AnyNIRSPEC) = read_key(data, "OBJECT")
Echelle.get_itime(data::AnyNIRSPEC) = read_key(data, "TRUITIME")
Echelle.get_image_num(data::AnyNIRSPEC) = parse(Int, split(split(basename(data), '_')[2], '.')[1])
Echelle.get_exposure_start_time(data::AnyNIRSPEC) = read_key(data, "MJD") + 2400000.5
Echelle.get_airmass(data::AnyNIRSPEC) = read_key(data, "AIRMASS")
Echelle.get_sky_coord(data::AnyNIRSPEC) = ICRSCoords(hms2rad(read_key(data, "RA")), dms2rad(read_key(data, "DEC")))
Echelle.get_timestamp(data::AnyNIRSPEC) = "20" * basename(data)[6:11]
Echelle.get_utdate(data::AnyNIRSPEC) = "20" * basename(data)[6:11]

# Orders for different modes
function get_orders(data::AnyNIRSPEC)
    m = get_mode(data)
    return ECHELLE_ORDERS[m]
end

# Orders for different modes
function get_mode(data::AnyNIRSPEC)
    echlpos = read_key(data, "ECHLPOS")
    disppos = read_key(data, "DISPPOS")
    return DISP_ANGLES_MODES[(echlpos, disppos)]
end

# Reading images
function Base.read(data::NIRSPECL0; hdu::Int=1, trans=false)
    image = FITS(data.filename) do f
        Float64.(read(f[hdu]))
    end
    if trans
        image .= collect(transpose(image))
    end
    return image
end

# mode-> angles
function get_angles(mode::String)
    angles = nothing
    for key in keys(DISP_ANGLES_MODES)
        if DISP_ANGLES_MODES[key] == mode
            angles = key
            break
        end
    end
    return angles
end