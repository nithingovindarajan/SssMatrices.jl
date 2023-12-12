function diagonalfill!(Adense, diagonal, block, no_blocks)
    for i = 1:no_blocks
        Adense[block[i], block[i]] = diagonal.D[i]
    end
end

function triangularfill!(Adense, triang, block, no_blocks)
    for j = 1:no_blocks
        temp = triang.inp[j]'
        for i = j+1:no_blocks
            Adense[block[i], block[j]] = triang.out[i] * temp
            temp = triang.trans[i] * temp
        end
    end
end

function Base.:Matrix(A::SSS)
    Adense = zeros(eltype(A), A.N, A.N)
    #fill up diagonal part
    diagonalfill!(Adense, A.diagonal, A.block, A.no_blocks)
    # fill up lower triangular part
    triangularfill!(Adense, A.lower, A.block, A.no_blocks)
    # fill up upper triangular part
    triangularfill!(Adense', A.upper, A.block, A.no_blocks)
    return Adense
end

function Base.:Matrix(A::Adjoint{Scalar,SSS{Scalar}}) where {Scalar<:Number}
    return Matrix(Matrix(A')')
end