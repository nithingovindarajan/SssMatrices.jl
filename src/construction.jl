function extract_diagonalpart(A, off, no_blocks)

    D = Matrix[]
    for i = 1:no_blocks
        push!(D, A[off[i]+1:off[i+1], off[i]+1:off[i+1]])
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
    Hbar = zeros(off[end] - off[2], 0)

    for i = 1:no_blocks-1

        Htilde = [Hbar A[off[i+1]+1:end, off[i]+1:off[i+1]]]

        # compute low rank approximation of Htilde
        X, Y = lowrankapprox(Htilde, threshold)     # Note there is no gaurantee that X and Y have entries in Scalar: potential source of bugs


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


function SSS_generators(A, n; threshold=1E-14)


    # assertions
    @assert size(A, 1) == size(A, 2) "Matrix is not square"
    @assert sum(n) == size(A, 1) "Matrix partition is not valid"

    # some definitions
    no_blocks = length(n)
    off = [0; cumsum(n)]

    # Define
    D = extract_diagonalpart(A, off, no_blocks)
    Q, R, P = extract_triangularpart(A, n, off, no_blocks, threshold)
    U, W, V = extract_triangularpart(A', n, off, no_blocks, threshold)

    return (D, Q, R, P, U, W, V, no_blocks, off)

end



function SSS{T}(A::AbstractMatrix, n::Vector{<:Integer}; threshold::Float64=1E-14) where {T<:Number}
    D, Q, R, P, U, W, V, no_blocks, off = SSS_generators(A, n; threshold=threshold)
    return SSS{T}(D, Q, R, P, U, W, V)
end





