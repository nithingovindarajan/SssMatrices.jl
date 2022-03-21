function diagonal_multiply!(B, diagonal::DiagonalPart, X, off::Vector{<:Integer}, no_blocks::Integer)

    for i = 1:no_blocks
        B[off[i]+1:off[i+1], :] += diagonal.D[i] * X[off[i]+1:off[i+1], :]
    end

end

function forward_iterate!(B, triang::TriangularPart, X, off::Vector{<:Integer}, no_blocks::Integer)

    # forward flow
    H = triang.inp[1]' * X[off[1]+1:off[2], :]
    for i = 2:no_blocks-1
        B[off[i]+1:off[i+1], :] += triang.out[i] * H
        H = triang.trans[i] * H + triang.inp[i]' * X[off[i]+1:off[i+1], :]
    end
    B[off[no_blocks]+1:off[no_blocks+1], :] += triang.out[no_blocks] * H

end


function backward_iterate!(B, triang::TriangularPart, X, off::Vector{<:Integer}, no_blocks::Integer)


    # backward flow
    G = triang.out[no_blocks]' * X[off[no_blocks]+1:off[no_blocks+1], :]
    for i = no_blocks-1:-1:2
        B[off[i]+1:off[i+1], :] += triang.inp[i] * G
        G = triang.trans[i]' * G + triang.out[i]' * X[off[i]+1:off[i+1], :]
    end
    B[off[1]+1:off[2], :] += triang.inp[1] * G


end


function SSSmultiply!(B, A::SSS, X)

    # multiply diagonal part
    diagonal_multiply!(B, A.diagonal, X, A.off, A.no_blocks)
    # multiply lower triangular part
    forward_iterate!(B, A.lower, X, A.off, A.no_blocks)
    # multiply upper triangular part
    backward_iterate!(B, A.upper, X, A.off, A.no_blocks)

end

function Base.:*(A::SSS, x::AbstractVector)


    # assertions
    @assert size(A, 2) == size(x, 1) "Matrix vector dimensions not consistent"

    b = zeros(size(A, 1))

    SSSmultiply!(b, A, x)

    return b

end




function Base.:*(A::SSS, X::AbstractMatrix)


    # assertions
    @assert size(A, 2) == size(X, 1) "Matrix vector dimensions not consistent"

    B = zeros(size(A, 1), size(X, 2))

    SSSmultiply!(B, A, X)

    return B

end
