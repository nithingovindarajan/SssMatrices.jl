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


    D, Q, R, P, U, W, V, block, no_blocks = SSS_generators_Cauchy(n; K=K, threshold=threshold)


    D = [sumofrowandcolumnscalings(D[l], view(Rtilde, block[l], :), view(Stilde, block[l], :)) for l ∈ 1:no_blocks]

    Q = [hcat(Tuple(Q[l] .* view(Stilde, block[l], i) for i ∈ 1:r)...) for l ∈ 1:no_blocks]
    R = [BlockDiagonal([R[l] for i ∈ 1:r]) for l ∈ 1:no_blocks]
    P = [hcat(Tuple(view(Rtilde, block[l], i) .* P[l] for i ∈ 1:r)...) for l ∈ 1:no_blocks]

    U = [hcat(Tuple(U[l] .* view(Rtilde, block[l], i) for i ∈ 1:r)...) for l ∈ 1:no_blocks]
    W = [BlockDiagonal([W[l] for i ∈ 1:r]) for l ∈ 1:no_blocks]
    V = [hcat(Tuple(view(Stilde, block[l], i) .* V[l] for i ∈ 1:r)...) for l ∈ 1:no_blocks]




    return SSS{ComplexF64}(D, Q, R, P, U, W, V)


end


#------------------------------------------------------------#
export SquareToeplitz


struct SquareToeplitz{Scalar<:Number} <: AbstractMatrix{Scalar}

    coeffs::Vector{Scalar}
    n::Int
    SSSform::SSS

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


