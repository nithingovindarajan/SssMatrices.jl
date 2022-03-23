module SSSmatrices

# dependent packages
using LinearAlgebra
using Random
using SparseArrays
using Tullio

export greet
greet() = print("Hello! welcome to the SSSmatrices package")


# include relevant files
include("SSStypes.jl")
include("sizes.jl")
include("matrixvectormultiplication.jl")
include("getindex.jl")
include("matrix.jl")
include("lra.jl")
include("construction.jl")
#include("addition.jl")
#include("matrixmatrixmultiplication.jl")
#include("LU.jl")
#include("ULV.jl")
#include("QR.jl")
#include("solve.jl")
include("CauchySSS.jl")
include("LDRmatrices.jl")
#include("casestudies.jl")











end # module
