function diagonal_multiply!(B, diagonal, X, block, no_blocks)

    for i = 1:no_blocks
        B[block[i], :] += diagonal.D[i] * X[block[i], :]
    end

end

function diagonal_multiply_adjoint!(B, diagonal, X, block, no_blocks)

    for i = 1:no_blocks
        B[block[i], :] += diagonal.D[i]' * X[block[i], :]
    end

end


function forward_iterate!(B, triang, X, block, no_blocks)

    # forward flow
    H = triang.inp[1]' * X[block[1], :]
    for i = 2:no_blocks-1
        B[block[i], :] += triang.out[i] * H
        H = triang.trans[i] * H + triang.inp[i]' * X[block[i], :]
    end
    B[block[no_blocks], :] += triang.out[no_blocks] * H

end


function backward_iterate!(B, triang, X, block, no_blocks)


    # backward flow
    G = triang.out[no_blocks]' * X[block[no_blocks], :]
    for i = no_blocks-1:-1:2
        B[block[i], :] += triang.inp[i] * G
        G = triang.trans[i]' * G + triang.out[i]' * X[block[i], :]
    end
    B[block[1], :] += triang.inp[1] * G


end



function Base.:*(A::SSS, x::AbstractVector)


    # assertions
    @assert size(A, 2) == size(x, 1) "Matrix vector dimensions not consistent"

    b = zeros(size(A, 1))

    # multiply diagonal part
    diagonal_multiply!(b, A.diagonal, x, A.block, A.no_blocks)
    # multiply lower triangular part
    forward_iterate!(b, A.lower, x, A.block, A.no_blocks)
    # multiply upper triangular part
    backward_iterate!(b, A.upper, x, A.block, A.no_blocks)

    return b

end


function Base.:*(A::SSS, X::AbstractMatrix)

    # assertions
    @assert size(A, 2) == size(X, 1) "Matrix vector dimensions not consistent"

    B = zeros(size(A, 1), size(X, 2))

    # multiply diagonal part
    diagonal_multiply!(B, A.diagonal, X, A.block, A.no_blocks)
    # multiply lower triangular part
    forward_iterate!(B, A.lower, X, A.block, A.no_blocks)
    # multiply upper triangular part
    backward_iterate!(B, A.upper, X, A.block, A.no_blocks)

    return B

end


function Base.:*(A::Adjoint{Scalar,SSS{Scalar}}, X::AbstractMatrix) where {Scalar<:Number}


    # assertions
    @assert size(A, 2) == size(X, 1) "Matrix vector dimensions not consistent"

    B = zeros(size(A, 1), size(X, 2))

    # multiply diagonal part
    diagonal_multiply_adjoint!(B, A.parent.diagonal, X, A.parent.block, A.parent.no_blocks)
    # multiply lower triangular part
    forward_iterate!(B, A.parent.upper, X, A.parent.block, A.parent.no_blocks)
    # multiply upper triangular part
    backward_iterate!(B, A.parent.lower, X, A.parent.block, A.parent.no_blocks)

    return B

end


function Base.:*(X::AbstractMatrix, A::Union{SSS,Adjoint{Scalar,SSS{Scalar}}}) where {Scalar<:Number}

    return (A' * X')'

end

