function SSS(A::Matrix, n::Vector{Int64}, threshold::Float64 = 1E-10)
    # assertions
    @assert size(A, 1) == size(A, 2) "Matrix is not square"
    @assert sum(n) == size(A, 1) "Matrix partition is not valid"


    # some definitions
    N = length(n)
    off = [0; cumsum(n)]

    # Define diagonal matrices
    Di = Matrix{Float64}[]
    for i = 1:N
        push!(Di, A[off[i]+1:off[i]+n[i], off[i]+1:off[i]+n[i]])
    end

    ### upper hankel blocks ###

    Hpi = Int64[]
    Ui = Matrix{Float64}[]
    Wi = Matrix{Float64}[]
    Vi = Matrix{Float64}[]

    # initialize
    Z = zeros(0, sum(n[2:end]))
    Un = zeros(0, 0)
    r = 0
    Leftprev = []
    for i = 1:N-1


        Z = vcat(Z, A[off[i]+1:off[i+1], off[i+1]+1:end])
        U, sigma, V, p = lowrankapprox(Z, threshold)
        Un = [Un * U[1:r, :]
            U[r+1:end, :]]
        r = p

        Left = Un * Diagonal(sigma)

        # add Hpi
        push!(Hpi, p)
        # add Ui
        push!(Ui, Left[off[i]+1:end, :])

        # add Vi
        push!(Vi, V[1:n[i+1], :])
        if i != 1
            # add Wi
            push!(Wi, Leftprev \ Left[1:off[i], :])
        end

        Z = Diagonal(sigma) * transpose(V[n[i+1]+1:end, :])
        Leftprev = deepcopy(Left)

    end

    ### lower hankel blocks ###

    Gpi = Int64[]
    Pi = Matrix{Float64}[]
    Ri = Matrix{Float64}[]
    Qi = Matrix{Float64}[]

    # initialize
    Z = zeros(sum(n[2:end]), 0)
    Vn = zeros(0, 0)
    r = 0
    Rightprev = []
    for i = 1:N-1

        Z = hcat(Z, A[off[i+1]+1:end, off[i]+1:off[i+1]])
        U, sigma, V, p = lowrankapprox(Z, threshold)
        Vn = [Vn * V[1:r, :]
            V[r+1:end, :]]
        r = p

        Right = Vn * Diagonal(sigma)


        # add Gpi
        push!(Gpi, p)
        # add Pi
        push!(Pi, U[1:n[i+1], :])
        # add Qi
        push!(Qi, Right[off[i]+1:end, :])
        if i != 1
            # add Ri
            push!(Ri, transpose(Rightprev \ Right[1:off[i], :]))
        end

        Z = U[n[i+1]+1:end, :] * Diagonal(sigma)
        Rightprev = deepcopy(Right)
    end


    #construct SSS matrix
    A_SSS = SSS{Float64}(N, n, Gpi, Hpi, Di, Ui, Wi, Vi, Pi, Ri, Qi)

    return A_SSS
end
