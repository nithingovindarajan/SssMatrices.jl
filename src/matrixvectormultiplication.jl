function diagonal_multiply!(b::Vector, diagonal::DiagonalPart, x::Vector, off::Vector{Integer}, N::Integer)

    for i = 1:N
        b[off[i]+1:block_i[i+1]] += diagonal[i] * x[off[i]+1:block_i[i+1]]
    end

end

function forward_iterate!(b::Vector, triang::TriangularPart, x::Vector, off::Vector{Integer}, N::Integer)

    # forward flow
    h = triang.inp[1]' * x[off[1]+1:off[2]]
    for i = 2:N-1
        b[off[i]+1:off[i+1]] += triang.out[i] * h
        h = triang.trans[i] * h + triang.inp[i]' * x[off[i]+1:off[i+1]]
    end
    b[off[N]+1:off[N+1]] += triang.out[N] * h

end


function backward_iterate!(b::Vector, triang::TriangularPart, x::Vector, off::Vector{Integer}, N::Integer)


    # backward flow
    g = triang.out[N]' * x[off[1]+1:off[2]]
    for i = N-1:-1:2
        b[off[i]+1:off[i+1]] += triang.inp[i] * g
        g = triang.trans[i]' * g + triang.inp[i] * xi[off[i]+1:off[i+1]]
    end
    b[off[1]+1:off[2]] += A.inp[1] * h


end

function Base.:*(A::SSS, x::Vector)
    # assertions
    @assert A.N == length(x) "Matrix vector dimensions not consistent"

    b = zeros(A.N)

    # multiply diagonal part
    diagonal_multiply!(b, A.diagonal, x, A.off, A.N)
    # multiply lower triangular part
    forward_iterate!(b, A.lower, x, A.off, A.N)
    # multiply upper triangular part
    backward_iterate!(b, A.upper, x, A.off, A.N)

    return b
end