function partitionVector(b::Vector, n::Vector{Int64})
    @assert sum(n) == length(b) "Matrix vector dimensions not consistent"
    off = [0; cumsum(n)]
    bi = Vector[b[off[i]+1:off[i]+n[i]] for i = 1:length(n)]

    return bi
end

function lowrankapprox(B::AbstractArray, threshold::Float64)

    # compute SVD
    U, sigma, V = svd(B)

    # rank
    p = findlast(x -> x > threshold, sigma)
    if p == nothing
        p = 0
    end

    #truncate
    U = U[:, 1:p]
    sigma = sigma[1:p]
    V = V[:, 1:p]

    return U, sigma, V, p
end

function eye(n::Integer)

    return Array(I, n, n)

end

