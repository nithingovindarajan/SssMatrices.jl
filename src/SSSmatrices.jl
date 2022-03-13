module SSSmatrices

# dependent packages
using LinearAlgebra
using Random
using SparseArrays

export greet
greet() = print("Hello! welcome to the SSSmatrices package")


# include relevant files
include("utilities.jl")
include("SSStypes.jl")
include("construction.jl")
include("addition.jl")
include("matrixvectormultiplication.jl")
include("matrixmatrixmultiplication.jl")
include("solve.jl")
include("LU.jl")
include("ULV.jl")
include("QR.jl")
include("casestudies.jl")











end # module
