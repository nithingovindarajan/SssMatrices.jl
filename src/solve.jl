
function Base.:\(A::SSS, b::Vector)


    N_sparse = sum(A.Gpi) + sum(A.Hpi) + sum(A.n)

    Blocksizes = Int64[]
    vecsizes = Array{Int64}(undef, A.N, 3)
    for i = 1:A.N
        if i == 1
            vecsizes[i, :] = [A.n[i] 0 A.Gpi[i]]
            append!(Blocksizes, A.n[i] + A.Gpi[i])
        elseif i == A.N
            vecsizes[i, :] = [A.n[i] A.Hpi[i-1] 0]
            append!(Blocksizes, A.n[i] + A.Hpi[i-1])
        else
            vecsizes[i, :] = [A.n[i] A.Hpi[i-1] A.Gpi[i]]
            append!(Blocksizes, A.n[i] + A.Gpi[i] + A.Hpi[i-1])
        end
    end

    Indices = pushfirst!(cumsum(Blocksizes), 0)
    pop!(Indices)

    #### ----------------------------------- ####
    #### ----  construct sparse vector  ---- ####
    #### ----------------------------------- ####


    I = Int64[]
    for i = 1:A.N
        append!(I, collect((Indices[i]+1):(Indices[i]+A.n[i])))
    end

    b_sparse = sparsevec(I, b, N_sparse)



    #### ----------------------------------- ####
    #### ----  construct sparse matrix  ---- ####
    #### ----------------------------------- ####

    I = Int64[]
    J = Int64[]
    V = Float64[]
    # add the A's
    for i = 1:A.N

        if i == 1
            Ak = [A.Di[i] zeros(vecsizes[i, 1], vecsizes[i, 2]) zeros(vecsizes[i, 1], vecsizes[i, 3])
                zeros(vecsizes[i, 2], vecsizes[i, 1]) -eye(vecsizes[i, 2]) zeros(vecsizes[i, 2], vecsizes[i, 3])
                transpose(A.Qi[i]) zeros(vecsizes[i, 3], vecsizes[i, 2]) -eye(vecsizes[i, 3])]
            Itemp, Jtemp, Vtemp = findnz(sparse(Ak))
        elseif i == A.N
            Ak = [A.Di[A.N] zeros(vecsizes[i, 1], vecsizes[i, 2]) zeros(vecsizes[i, 1], vecsizes[i, 3])
                transpose(A.Vi[A.N-1]) -eye(vecsizes[i, 2]) zeros(vecsizes[i, 2], vecsizes[i, 3])
                zeros(vecsizes[i, 3], vecsizes[i, 1]) zeros(vecsizes[i, 3], vecsizes[i, 2]) -eye(vecsizes[i, 3])]
            Itemp, Jtemp, Vtemp = findnz(sparse(Ak))
        else
            Ak = [A.Di[i] zeros(vecsizes[i, 1], vecsizes[i, 2]) zeros(vecsizes[i, 1], vecsizes[i, 3])
                transpose(A.Vi[i-1]) -eye(vecsizes[i, 2]) zeros(vecsizes[i, 2], vecsizes[i, 3])
                transpose(A.Qi[i]) zeros(vecsizes[i, 3], vecsizes[i, 2]) -eye(vecsizes[i, 3])]
            Itemp, Jtemp, Vtemp = findnz(sparse(Ak))
        end

        append!(I, Indices[i] .+ Itemp)
        append!(J, Indices[i] .+ Jtemp)
        append!(V, Vtemp)


    end


    # add the B's
    for i = 2:A.N
        if i == A.N
            Bk = [zeros(vecsizes[i, 1], vecsizes[i-1, 1]) zeros(vecsizes[i, 1], vecsizes[i-1, 2]) A.Pi[i-1]
                zeros(vecsizes[i, 2], vecsizes[i-1, 1]) zeros(vecsizes[i, 2], vecsizes[i-1, 2]) zeros(vecsizes[i, 2], vecsizes[i-1, 3])
                zeros(vecsizes[i, 3], vecsizes[i-1, 1]) zeros(vecsizes[i, 3], vecsizes[i-1, 2]) zeros(vecsizes[i, 3], vecsizes[i-1, 3])]
            Itemp, Jtemp, Vtemp = findnz(sparse(Bk))
        else
            Bk = [zeros(vecsizes[i, 1], vecsizes[i-1, 1]) zeros(vecsizes[i, 1], vecsizes[i-1, 2]) A.Pi[i-1]
                zeros(vecsizes[i, 2], vecsizes[i-1, 1]) zeros(vecsizes[i, 2], vecsizes[i-1, 2]) zeros(vecsizes[i, 2], vecsizes[i-1, 3])
                zeros(vecsizes[i, 3], vecsizes[i-1, 1]) zeros(vecsizes[i, 3], vecsizes[i-1, 2]) A.Ri[i-1]]
            Itemp, Jtemp, Vtemp = findnz(sparse(Bk))
        end

        append!(I, Indices[i] .+ Itemp)
        append!(J, Indices[i-1] .+ Jtemp)
        append!(V, Vtemp)


    end

    # add the C's
    for i = 2:A.N
        if i == 2
            Ck = [zeros(vecsizes[i-1, 1], vecsizes[i, 1]) A.Ui[i-1] zeros(vecsizes[i-1, 1], vecsizes[i, 3])
                zeros(vecsizes[i-1, 2], vecsizes[i, 1]) zeros(vecsizes[i-1, 2], vecsizes[i, 2]) zeros(vecsizes[i-1, 2], vecsizes[i, 3])
                zeros(vecsizes[i-1, 3], vecsizes[i, 1]) zeros(vecsizes[i-1, 3], vecsizes[i, 2]) zeros(vecsizes[i-1, 3], vecsizes[i, 3])]
            Itemp, Jtemp, Vtemp = findnz(sparse(Ck))
        else
            Ck = [zeros(vecsizes[i-1, 1], vecsizes[i, 1]) A.Ui[i-1] zeros(vecsizes[i-1, 1], vecsizes[i, 3])
                zeros(vecsizes[i-1, 2], vecsizes[i, 1]) A.Wi[i-2] zeros(vecsizes[i-1, 2], vecsizes[i, 3])
                zeros(vecsizes[i-1, 3], vecsizes[i, 1]) zeros(vecsizes[i-1, 3], vecsizes[i, 2]) zeros(vecsizes[i-1, 3], vecsizes[i, 3])]
            Itemp, Jtemp, Vtemp = findnz(sparse(Ck))
        end

        append!(I, Indices[i-1] .+ Itemp)
        append!(J, Indices[i] .+ Jtemp)
        append!(V, Vtemp)


    end

    A_sparse = sparse(I, J, V, N_sparse, N_sparse)

    #### ------------------------------- ####
    #### ----  solve sparse system  ---- ####
    #### ------------------------------- ####

    x_sparse = qr(A_sparse) \ Vector(b_sparse)

    x = Float64[]
    for i = 1:A.N
        append!(x, x_sparse[(Indices[i]+1):(Indices[i]+A.n[i])])
    end

    return x

end
