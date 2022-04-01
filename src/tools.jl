export offdiagonal_num_rank

function offdiagonal_num_rank(A, n; threshold=1E-12)
    @assert size(A, 1) == size(A, 2)
    # compute SVD
    U, sigma, V = svd(A[n+1:size(A, 1), 1:n])

    # rank
    p = findlast(x -> x > sigma[1] * threshold, sigma)
    if p == nothing
        p = 0
    end

    return p
end



export rank_lower_hankelblock
function rank_lower_hankelblock(A, n, k; threshold=1E-12)

    @assert size(A, 1) == size(A, 2)


    # some definitions
    no_blocks = length(n)
    off = [0; cumsum(n)]

    # compute SVD
    U, sigma, V = svd(A[off[k+1]+1:end, 1:off[k+1]])

    # rank
    p = findlast(x -> x > sigma[1] * threshold, sigma)
    if p == nothing
        p = 0
    end

    return p
end
