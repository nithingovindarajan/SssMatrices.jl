function Base.:*(A::SSS, x::Vector)
    # assertions
    @assert sum(A.n) == length(x) "Matrix vector dimensions not consistent"

    # break up vector
    xi = partitionVector(x, A.n)

    N = A.N

    ### SSS multiply ###
    # Diagonal terms
    bi = [A.Di[i] * xi[i] for i = 1:N]
    # backward flow (upper triangular part)
    g = transpose(A.Vi[N-1]) * xi[N]
    bi[N-1] = bi[N-1] + A.Ui[N-1] * g
    for i = N-2:-1:1
        g = A.Wi[i] * g + transpose(A.Vi[i]) * xi[i+1]
        bi[i] = bi[i] + A.Ui[i] * g
    end
    # forward flow (lower triangular part)
    h = transpose(A.Qi[1]) * xi[1]
    bi[2] = bi[2] + A.Pi[1] * h
    for i = 2:1:N-1
        h = A.Ri[i-1] * h + transpose(A.Qi[i]) * xi[i]
        bi[i+1] = bi[i+1] + A.Pi[i] * h
    end

    #concat back into vector
    b = foldr(vcat, bi)

    return b
end