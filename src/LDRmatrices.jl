
Dminone(n::Int) = Diagonal([(-1.0 + 0im)^(-(k - 1) / n) for k = 1:n])
OMEGA(n::Int) = [exp(π * 2im * (k - 1) / n) for k = 1:n]
LAMBDA(n::Int) = exp(π * 1im / n) * OMEGA(n::Int)

#------------------------------------------------------------#

export CauchyLike


struct CauchyLike{Scalar<:Number} <: AbstractMatrix{Scalar}
    n::Integer
    r::Integer
    omega::Vector{Scalar}
    lambda::Vector{Scalar}
    R::Matrix{Scalar}
    S::Matrix{Scalar}

    function CauchyLike{Scalar}(omega, lambda, R, S) where {Scalar<:Number}

        if length(omega) == length(lambda) == size(R, 1) == size(S, 1) && size(R, 2) == size(S, 2)

            n = length(omega)
            r = size(R, 2)
            new{Scalar}(n, r, convert(Vector{Scalar}, omega), convert(Vector{Scalar}, lambda),
                convert(Matrix{Scalar}, R), convert(Matrix{Scalar}, S))

        else
            DomainError("Provided inputs have incompatible dimensions")
        end


    end

end
CauchyLike{Scalar}(omega, lambda) where {Scalar<:Number} = CauchyLike{Scalar}(omega, lambda, ones(Scalar, length(omega), 1), ones(Scalar, length(omega), 1))

Base.:size(A::CauchyLike) = (A.n, A.n)

function Base.:getindex(A::CauchyLike, i::Int, j::Int)

    if 1 <= i <= A.n && 1 <= j <= A.n
        return dot(A.S[j, :], A.R[i, :]) / (A.omega[i] - A.lambda[j])
    else
        BoundsError()
    end

end

#------------------------------------------------------------#



# must return the SSS representation of the Cauchy matrix of the Cauchy trick
function SSS_CauchyLike(Rtilde, Stilde)

    n = size(Rtilde, 1)

    return CauchyLike{ComplexF64}(OMEGA(n), LAMBDA(n), Rtilde, Stilde) ## temporary


end


#------------------------------------------------------------#
export SquareToeplitz

# determines the corresponding chauchy matrix (in SSS form) of the Toeplitz matrix in SS
function ToeplitzSSSform(Rtilde, Stilde)

    n = size(Rtilde, 1)

    return CauchyLike{ComplexF64}(OMEGA(n), LAMBDA(n), Rtilde, Stilde) ## temporary


end

struct SquareToeplitz{Scalar<:Number} <: AbstractMatrix{Scalar}

    coeffs::Vector{Scalar}
    n::Int
    SSSform::CauchyLike{<:Complex}    # To be replaced with SSS{Scalar} !!!

    function SquareToeplitz(coeffs::Vector{T}, n::Integer) where {T<:Number}
        if length(coeffs) == 2 * n - 1


            R = [1 coeffs[n]; zeros(n - 1) coeffs[1:n-1]+coeffs[n+1:end]]
            S = [coeffs[end:-1:n+1]-coeffs[n-1:-1:1] zeros(n - 1); coeffs[n] 1]
            Rtilde = sqrt(n) * ifft(R, 1)
            Stilde = sqrt(n) * ifft((Dminone(n))' * S, 1)


            new{T}(coeffs, n, SSS_CauchyLike(Rtilde, Stilde))
        else
            DimensionMismatch()
        end
    end

end
function SquareToeplitz(coeff_lt::Vector, coeff_ut::Vector) where {Scalar<:Number}
    n = length(coeff_lt)
    coeffs = [reverse(coeff_ut); coeff_lt]
    return SquareToeplitz(coeffs, n)
end
Base.:size(A::SquareToeplitz) = (A.n, A.n)
function Base.:getindex(A::SquareToeplitz, i::Int, j::Int)

    if 1 <= i <= A.n && 1 <= j <= A.n
        return A.coeffs[i-j+A.n]
    else
        BoundsError()
    end

end


function Base.:\(A::SquareToeplitz, b::Vector)


    if length(b) != size(A, 2)
        DimensionMismatch()
    else

        btilde = sqrt(A.n) * ifft(b)
        xtilde = \(A.SSSform, btilde)
        x = Dminone(A.n) * (1 / sqrt(A.n)) * fft(xtilde)

        return x
    end

end
