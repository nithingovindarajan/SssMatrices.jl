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
