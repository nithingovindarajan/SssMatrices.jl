function diagonalfill!(A, diagonal::DiagonalPart, off::Vector{Integer}, N::Integer)

    for i = 1:N
        A[off[i]+1:block_i[i+1], off[i]+1:block_i[i+1]] = diagonal[i]
    end

end


function triangularfill!(A, triang::TriangularPart, off::Vector{Integer}, N::Integer)

    for j = 1:N

        temp = triang.inp[j]'

        for i = j+1:n

            Adense[off[i]+1:block_i[i+1], off[j]+1:block_i[j+1]] = triang.out[i] * temp
            temp = triang.trans[i] * temp

        end

    end

end


function Base.:Matrix{Scalar}(A::SSS) where {Scalar<:Number}


    Adense = zeros(Scalar, ntot, ntot)

    #fill up diagonal part
    diagonalfill!(Adense, A.diagonal, A.off, A.N)
    # fill up lower triangular part
    triangularfill!(Adense, A.lower, A.off, A.N)
    # fill up upper triangular part
    triangularfill!(Adense', A.upper, A.off, A.N)

    return Adense
end


