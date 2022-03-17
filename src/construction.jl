function extract_diagonalpart(A::Matrix, n::Vector{Int}, N::Int, off:::Vector{Int})

    D = Matrix{Float64}[]
    for i = 1:N
        push!(D, A[off[i]+1:off[i]+n[i], off[i]+1:off[i]+n[i]])
    end

    return D

end


function extract_triangularpart(A::Matrix, n::Vector{Int}, N::Int) where {Scalar<:Number}

    Scalar = eltype(A)
    Pi = Matrix{Scalar}[]
    Ri = Matrix{Scalar}[]
    Qi = Matrix{Scalar}[]

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

    return Q, R, P


end


function SSS(A::Matrix, n::Vector{Int64}, threshold::Float64 = 1E-10)


    # assertions
    @assert size(A, 1) == size(A, 2) "Matrix is not square"
    @assert sum(n) == size(A, 1) "Matrix partition is not valid"

    # some definitions
    N = length(n)
    off = [0; cumsum(n)]

    # Define
    D = extract_diagonalpart(A, n, N, off)
    Q, R, P = extract_triangularpart(A, n, N, off)
    U, W, V = extract_triangularpart(A', n, N, off)

    return SSS(D, Q, R, P, U, W, V)

end
