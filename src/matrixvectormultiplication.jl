function diagonal_multiply!(B, diagonal, X, off, no_blocks)

    for i = 1:no_blocks
        B[off[i]+1:off[i+1], :] += diagonal.D[i] * X[off[i]+1:off[i+1], :]
    end

end

function diagonal_multiply_adjoint!(B, diagonal, X, off, no_blocks)

    for i = 1:no_blocks
        B[off[i]+1:off[i+1], :] += diagonal.D[i]' * X[off[i]+1:off[i+1], :]
    end

end


function forward_iterate!(B, triang, X, off, no_blocks)

    # forward flow
    H = triang.inp[1]' * X[off[1]+1:off[2], :]
    for i = 2:no_blocks-1
        B[off[i]+1:off[i+1], :] += triang.out[i] * H
        H = triang.trans[i] * H + triang.inp[i]' * X[off[i]+1:off[i+1], :]
    end
    B[off[no_blocks]+1:off[no_blocks+1], :] += triang.out[no_blocks] * H

end


function backward_iterate!(B, triang, X, off, no_blocks)


    # backward flow
    G = triang.out[no_blocks]' * X[off[no_blocks]+1:off[no_blocks+1], :]
    for i = no_blocks-1:-1:2
        B[off[i]+1:off[i+1], :] += triang.inp[i] * G
        G = triang.trans[i]' * G + triang.out[i]' * X[off[i]+1:off[i+1], :]
    end
    B[off[1]+1:off[2], :] += triang.inp[1] * G


end



function Base.:*(A::SSS, x::AbstractVector)


    # assertions
    @assert size(A, 2) == size(x, 1) "Matrix vector dimensions not consistent"

    b = zeros(size(A, 1))

    # multiply diagonal part
    diagonal_multiply!(b, A.diagonal, x, A.off, A.no_blocks)
    # multiply lower triangular part
    forward_iterate!(b, A.lower, x, A.off, A.no_blocks)
    # multiply upper triangular part
    backward_iterate!(b, A.upper, x, A.off, A.no_blocks)

    return b

end


function Base.:*(A::SSS, X::AbstractMatrix)

    # assertions
    @assert size(A, 2) == size(X, 1) "Matrix vector dimensions not consistent"

    B = zeros(size(A, 1), size(X, 2))

    # multiply diagonal part
    diagonal_multiply!(B, A.diagonal, X, A.off, A.no_blocks)
    # multiply lower triangular part
    forward_iterate!(B, A.lower, X, A.off, A.no_blocks)
    # multiply upper triangular part
    backward_iterate!(B, A.upper, X, A.off, A.no_blocks)

    return B

end


function Base.:*(A::Adjoint{Scalar,SSS{Scalar}}, X::AbstractMatrix) where {Scalar<:Number}


    # assertions
    @assert size(A, 2) == size(X, 1) "Matrix vector dimensions not consistent"

    B = zeros(size(A, 1), size(X, 2))

    # multiply diagonal part
    diagonal_multiply_adjoint!(B, A.parent.diagonal, X, A.parent.off, A.parent.no_blocks)
    # multiply lower triangular part
    forward_iterate!(B, A.parent.upper, X, A.parent.off, A.parent.no_blocks)
    # multiply upper triangular part
    backward_iterate!(B, A.parent.lower, X, A.parent.off, A.parent.no_blocks)

    return B

end


function Base.:*(X::AbstractMatrix, A::Union{SSS,Adjoint{Scalar,SSS{Scalar}}}) where {Scalar<:Number}

    return (A' * X')'

end

