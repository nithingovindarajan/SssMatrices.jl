function diagonalfill!(Adense, diagonal::DiagonalPart, off::Vector{<:Integer}, no_blocks::Integer)

    for i = 1:no_blocks
        Adense[off[i]+1:off[i+1], off[i]+1:off[i+1]] = diagonal.D[i]
    end

end


function triangularfill!(Adense, triang::TriangularPart, off::Vector{<:Integer}, no_blocks::Integer)

    for j = 1:no_blocks

        temp = triang.inp[j]'

        for i = j+1:no_blocks

            Adense[off[i]+1:off[i+1], off[j]+1:off[j+1]] = triang.out[i] * temp
            temp = triang.trans[i] * temp

        end

    end

end


function Base.:Matrix{Scalar}(A::SSS) where {Scalar<:Number}


    Adense = zeros(Scalar, A.N, A.N)

    #fill up diagonal part
    diagonalfill!(Adense, A.diagonal, A.off, A.no_blocks)
    # fill up lower triangular part
    triangularfill!(Adense, A.lower, A.off, A.no_blocks)
    # fill up upper triangular part
    triangularfill!(Adense', A.upper, A.off, A.no_blocks)

    return Adense
end


