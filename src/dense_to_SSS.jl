function lowrankapprox(B::AbstractMatrix, threshold::Float64)
    # compute SVD
    U, sigma, V = svd(B)
    # rank
    p = findlast(x -> x > sigma[1] * threshold, sigma)
    if p == nothing
        p = 0
    end
    #truncate
    X = U[:, 1:p]
    Y = V[:, 1:p] * Diagonal(sigma[1:p])
    return X, Y, p
end

function extract_diagonalpart(A, block, no_blocks)
    D = Matrix[]
    for i = 1:no_blocks
        push!(D, A[block[i], block[i]])
    end
    return D
end

function extract_triangularpart(A, n, off, no_blocks, threshold)
    inp = Matrix[]
    trans = Matrix[]
    out = Matrix[]
    # add out_(1)
    push!(out, zeros(n[1], 0))
    # initialize Hbar
    Hbar = zeros(size(A, 2) - n[1], 0)
    for i = 1:no_blocks-1
        Htilde = [Hbar A[off[i+1]+1:end, off[i]+1:off[i+1]]]
        # compute low rank approximation of Htilde
        X, Y, p = lowrankapprox(Htilde, threshold)
        # add inp_(i)
        push!(inp, Y[size(Hbar, 2)+1:end, :])
        # add trans_(i)
        push!(trans, Y[1:size(Hbar, 2), :]')
        # add out_(i+1)
        push!(out, X[1:n[i+1], :])

        # determine next Hbar
        Hbar = X[n[i+1]+1:end, :]
    end
    # add inp_(N)
    push!(inp, zeros(n[no_blocks], 0))
    # add trans_(N)
    push!(trans, zeros(size(Hbar, 2), 0)')

    return inp, trans, out
end

function SSS_generators(A, n; threshold = 1E-14)
    # assertions
    @assert size(A, 1) == size(A, 2) "Matrix is not square"
    @assert sum(n) == size(A, 1) "Matrix partition is not valid"
    # some definitions
    no_blocks = length(n)
    off = [0; cumsum(n)]
    block = [off[i]+1:off[i+1] for i = 1:no_blocks]
    # Define
    D = extract_diagonalpart(A, block, no_blocks)
    Q, R, P = extract_triangularpart(A, n, off, no_blocks, threshold)
    U, W, V = extract_triangularpart(A', n, off, no_blocks, threshold)
    return (D, Q, R, P, U, W, V, block, no_blocks)
end

function SSS{T}(
    A::AbstractMatrix,
    n::Vector{<:Integer};
    threshold::Float64 = 1E-14,
) where {T<:Number}
    D, Q, R, P, U, W, V, _, _ = SSS_generators(A, n; threshold = threshold)
    return SSS{T}(D, Q, R, P, U, W, V)
end
