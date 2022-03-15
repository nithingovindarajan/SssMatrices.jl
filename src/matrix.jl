function Base.:Matrix{Scalar}(A::SSS) where {Scalar<:Number}
    # assertions

    ntot = sum(A.n)
    Adense = zeros(Scalar, ntot, ntot)
    off = [0; cumsum(A.n)]

    #fill up diagonal part
    for i = 1:A.N
        Adense[off[i]+1:off[i]+A.n[i], off[i]+1:off[i]+A.n[i]] = A.Di[i]
    end

    # fill up lower triangular part
    for j = 1:(A.N-1)
        temp = transpose(A.Qi[j])
        for i = j+1:A.N
            Adense[off[i]+1:off[i]+A.n[i], off[j]+1:off[j]+A.n[j]] = A.Pi[i-1] * temp
            if i != A.N
                temp = A.Ri[i-1] * temp
            end
        end
    end

    # fill up upper triangular part
    for j = A.N:-1:2
        temp = transpose(A.Vi[j-1])
        for i = j-1:-1:1
            Adense[off[i]+1:off[i]+A.n[i], off[j]+1:off[j]+A.n[j]] = A.Ui[i] * temp
            if i != 1
                temp = A.Wi[i-1] * temp
            end
        end
    end

    return Adense
end











