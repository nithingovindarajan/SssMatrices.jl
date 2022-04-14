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

export determine_blocksizes_cauchy

function determine_blocksizes_cauchy(n; K=1.0)

    size_block = Int(ceil((K * log2(n))))
    no_blocks, remainder = divrem(n, size_block)
    if remainder != 0
        return [fill(size_block, no_blocks); remainder]
    else
        return fill(size_block, no_blocks)
    end
end


# must return the SSS representation of the Cauchy matrix of the Cauchy trick. Function performance to be improved later...
function SSS_Cauchy(n; K=1.0, threshold=1E-14)

    return SSS{ComplexF64}(FourierCauchy(n), determine_blocksizes_cauchy(n; K=K); threshold=threshold)

end


SSS_generators_Cauchy(n; K=1.0, threshold=1E-14) = SSS_generators(FourierCauchy(n), determine_blocksizes_cauchy(n; K=K); threshold=threshold)
