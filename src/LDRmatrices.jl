
Dminone(n::Int) = Diagonal([(-1.0 + 0im)^(-(k - 1) / n) for k = 1:n])
OMEGA(n::Int) = [exp(π * 2im * (k - 1) / n) for k = 1:n]
LAMBDA(n::Int) = [exp(π * (2 * (k - 1) + 1) * im / n) for k = 1:n]


#------------------------------------------------------------#

export FourierCauchy

struct FourierCauchy{Scalar<:Complex{<:AbstractFloat}} <: AbstractMatrix{Scalar}
    n::Integer
    omega::Vector{Scalar}
    lambda::Vector{Scalar}


    function FourierCauchy{Scalar}(n::Integer) where {Scalar<:Complex{<:AbstractFloat}}

        if n > 0

            omega = convert(Vector{Scalar}, OMEGA(n))
            lambda = convert(Vector{Scalar}, LAMBDA(n))

            new{Scalar}(n, omega, lambda)
        else
            DomainError("please provide a positive integer")
        end


    end

end
FourierCauchy(n::Integer) = FourierCauchy{ComplexF64}(n::Integer)
Base.:size(A::FourierCauchy) = (A.n, A.n)
Base.:getindex(A::FourierCauchy, i::Int, j::Int) = 1 / (A.omega[i] - A.lambda[j])


#------------------------------------------------------------#

export SSS_Cauchy, determine_blocksizes

function determine_blocksizes(n, K)

    size_block = Int(ceil((K * log2(n))))
    no_blocks, remainder = divrem(n, size_block)

    return [fill(size_block, no_blocks); remainder]
end


# must return the SSS representation of the Cauchy matrix of the Cauchy trick. Function performance to be improved later...
function SSS_Cauchy(n; K=1.0, threshold=1E-14)

    return SSS(FourierCauchy(n), determine_blocksizes(n, K); threshold=threshold)

end


SSS_generators_Cauchy(n; K=1.0, threshold=1E-14) = SSS_generators(FourierCauchy(n), determine_blocksizes(n, K); threshold=threshold)

function sumofrowandcolumnscalings(D, Rtilde, Stilde)
    @tullio Dnew[i, j] := Rtilde[i, k] * D[i, j] * conj(Stilde[j, k])
    return Dnew
end

function SSS_CauchyLike(coeffs, n; K=1.0, threshold=1E-14)


    # determine Rtilde and Stilde
    r = 2
    Rmat = [1 coeffs[n]; zeros(n - 1) coeffs[1:n-1]+coeffs[n+1:end]]
    Smat = [coeffs[end:-1:n+1]-coeffs[n-1:-1:1] zeros(n - 1); coeffs[n] 1]
    Rtilde = sqrt(n) * ifft(Rmat, 1)
    Stilde = sqrt(n) * ifft((Dminone(n))' * Smat, 1)


    D, Q, R, P, U, W, V, no_blocks, off = SSS_generators_Cauchy(n; K=K, threshold=threshold)


    D = [sumofrowandcolumnscalings(D[l], view(Rtilde, off[l]+1:off[l+1], :), view(Stilde, off[l]+1:off[l+1], :)) for l ∈ 1:no_blocks]

    Q = [hcat(Tuple(Q[l] .* view(Stilde, off[l]+1:off[l+1], i) for i ∈ 1:r)...) for l ∈ 1:no_blocks]
    R = [BlockDiagonal([R[l] for i ∈ 1:r]) for l ∈ 1:no_blocks]
    P = [hcat(Tuple(view(Rtilde, off[l]+1:off[l+1], i) .* P[l] for i ∈ 1:r)...) for l ∈ 1:no_blocks]

    U = [hcat(Tuple(U[l] .* view(Rtilde, off[l]+1:off[l+1], i) for i ∈ 1:r)...) for l ∈ 1:no_blocks]
    W = [BlockDiagonal([W[l] for i ∈ 1:r]) for l ∈ 1:no_blocks]
    V = [hcat(Tuple(view(Stilde, off[l]+1:off[l+1], i) .* V[l] for i ∈ 1:r)...) for l ∈ 1:no_blocks]




    return SSS(D, Q, R, P, U, W, V)


end


#------------------------------------------------------------#
export SquareToeplitz


struct SquareToeplitz{Scalar<:Number} <: AbstractMatrix{Scalar}

    coeffs::Vector{Scalar}
    n::Int
    SSSform::AbstractMatrix   # To be replaced with SSS{Scalar} !!!

    function SquareToeplitz(coeffs::Vector{T}, n::Integer) where {T<:Number}
        if length(coeffs) == 2 * n - 1
            new{T}(coeffs, n, SSS_CauchyLike(coeffs, n))
        else
            DimensionMismatch()
        end
    end

end
function SquareToeplitz(coeff_lt::Vector, coeff_ut::Vector)
    n = length(coeff_lt)
    coeffs = [reverse(coeff_ut); coeff_lt]
    return SquareToeplitz(coeffs, n)
end
Base.:size(A::SquareToeplitz) = (A.n, A.n)
Base.:getindex(A::SquareToeplitz, i::Int, j::Int) = A.coeffs[i-j+A.n]



function Base.:\(A::SquareToeplitz, b::Vector)


    if length(b) != size(A, 2)
        DimensionMismatch()
    else

        btilde = sqrt(A.n) * ifft(b)
        xtilde = \(Matrix(A.SSSform), btilde)
        x = Dminone(A.n) * (1 / sqrt(A.n)) * fft(xtilde)

        return x
    end

end


