module SssMatrices

# dependent packages
using LinearAlgebra
using Random
using BlockDiagonals
using FFTW
using Tullio

# export 
export SSS
export random_SSS

include("SSStypes.jl")
include("addition.jl")
include("SSSmatvec.jl")
include("SSS_to_dense.jl")
include("dense_to_SSS.jl")
include("solve.jl")











end # module
