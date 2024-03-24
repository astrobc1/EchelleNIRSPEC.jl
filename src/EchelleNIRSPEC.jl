module EchelleNIRSPEC

using FITSIO, JLD2
using Infiltrator
using Polynomials
using OrderedCollections
using AstroAngles, SkyCoords
using NaNStatistics
using Glob
using Reexport
using Distributed

@reexport using Echelle
using EchelleReduce

# Exports
export NIRSPECL0, NIRSPECL1

# Basic info
const NAME = "NIRSPEC"
const OBSERVATORY = "Keck"
const TIMEZONE = -10

# Detector
const DETECTOR_GAIN = 1
const DETECTOR_READ_NOISE = 0
const DETECTOR_DARK_CURRENT = 0

const PATHSEP = Base.Filesystem.path_separator

const ECHELLE_ORDERS0 = Dict(
    "KCEN" => 1:4,
    "K1" => 1:9,
    "K2" => 1:9,
    "K3" => 1:6,
    "K4" => 1:5,
    "H" => 1:11,
    "HK" => 1:10,
    "J" => 1:8,
    "Y" => 1:6,
    "JH" => 1:1,
    "YJ" => 1:11
)

const ECHELLE_ORDERS = Dict(
    "Y" => 75:-1:68,
    "YJ" => 70:-1:60,
    "J" => 65:-1:56,
    # "JH" => 1:1, ??
    "H" => 53:-1:43,
    "HK" => 45:-1:35,
)

const DISP_ANGLES_MODES = Dict(
    (63.0, 35.76) => "KCEN",
    (62.22, 35.0) => "K1",
    (63.90, 35.17) => "K2",
    (62.15, 37.0) => "K3",
    (64.0, 37.2) => "K4",
    (63.0, 36.72) => "H",
    (63.0, 34.31) => "HK",
    (63.0, 34.08) => "J",
    (63.00, 34.95) => "Y",
    (63.00, 35.23) => "JH",
    (63.00, 36.48) => "YJ"
)

# Data types for nirspec -> L0 and L1 for all images and extract spectra, respectively.
abstract type AnyNIRSPEC{L} <: SpecData{Symbol(lowercase(NAME)), L} end
struct NIRSPECL0 <: AnyNIRSPEC{0}
    filename::String
end
struct NIRSPECL1 <: AnyNIRSPEC{1}
    filename::String
end


include("parsing.jl")


include("reduction.jl")


end
