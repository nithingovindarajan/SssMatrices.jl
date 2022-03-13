export DiagonalSSS, SSS


struct DiagonalSSS{Scalar<:Number} <: AbstractMatrix{Scalar}
    N::Int64
    n::Vector{Int64}
    Gpi::Vector{Int64}   # hankel block ranks lower triangular part
    Hpi::Vector{Int64}   # hankel block rank upper triangular part
    # main diagonal
    Di::Vector{Matrix{Scalar}}
    # upper triangular part
    Ui::Vector{Matrix{Scalar}}
    Wi::Vector{Matrix{Scalar}}
    Vi::Vector{Matrix{Scalar}}
    # lower triangular part
    Pi::Vector{Matrix{Scalar}}
    Ri::Vector{Matrix{Scalar}}
    Qi::Vector{Matrix{Scalar}}

    function DiagonalSSS{Scalar}(N, n, Gpi, Hpi, Di, Ui, Wi, Vi, Pi, Ri, Qi) where {Scalar<:Number}

        if N < 3
            error("N<3 not supported")
        end

        # check input dimensions
        if length(n) != N || length(Di) != N || length(Ui) != N - 1 ||
           length(Wi) != N - 2 || length(Vi) != N - 1 || length(Pi) != N - 1 ||
           length(Ri) != N - 2 || length(Qi) != N - 1 || length(Hpi) != N - 1 ||
           length(Gpi) != N - 1
            error("Input dimensions do not agree")
        end

        # check diagonal matrix dimensions
        for i = 1:N
            if size(Di[i], 1) != n[i] || size(Di[i], 2) != n[i]
                error("Diagonal matrices dimensions mismatch")
            end
        end

        # check upper triangular matrix dimensions
        for i = 1:N-1
            if size(Ui[i], 1) != n[i] || size(Ui[i], 2) != Hpi[i]
                error("Ui matrices dimensions mismatch")
            end
            if size(Vi[i], 1) != n[i+1] || size(Vi[i], 2) != Hpi[i]
                error("Vi matrices dimensions mismatch")
            end
        end
        for i = 1:N-2
            if size(Wi[i], 1) != Hpi[i] || size(Wi[i], 2) != Hpi[i+1]
                error("Wi matrices dimensions mismatch")
            end
        end

        # check lower triangular matrix dimensions
        for i = 1:N-1
            if size(Pi[i], 1) != n[i+1] || size(Pi[i], 2) != Gpi[i]
                error("Pi matrices dimensions mismatch")
            end
            if size(Qi[i], 1) != n[i] || size(Qi[i], 2) != Gpi[i]
                error("Qi matrices dimensions mismatch")
            end
        end
        for i = 1:N-2
            if size(Ri[i], 1) != Gpi[i+1] || size(Ri[i], 2) != Gpi[i]
                error("Ri translation operators dimensions mismatch")
            end
        end

        new{Scalar}(N, n, Gpi, Hpi,
            map(x -> convert(Matrix{Scalar}, x), Di),
            map(x -> convert(Matrix{Scalar}, x), Ui),
            map(x -> convert(Matrix{Scalar}, x), Wi),
            map(x -> convert(Matrix{Scalar}, x), Vi),
            map(x -> convert(Matrix{Scalar}, x), Pi),
            map(x -> convert(Matrix{Scalar}, x), Ri),
            map(x -> convert(Matrix{Scalar}, x), Qi))


    end

end
Base.:size(A::DiagonalSSS) = (sum(A.n), sum(A.n))


struct SSS{Scalar<:Number} <: AbstractMatrix{Scalar}
    N::Int64
    n::Vector{Int64}
    Gpi::Vector{Int64}   # hankel block ranks lower triangular part
    Hpi::Vector{Int64}   # hankel block rank upper triangular part
    # main diagonal
    Di::Vector{Matrix{Scalar}}
    # upper triangular part
    Ui::Vector{Matrix{Scalar}}
    Wi::Vector{Matrix{Scalar}}
    Vi::Vector{Matrix{Scalar}}
    # lower triangular part
    Pi::Vector{Matrix{Scalar}}
    Ri::Vector{Matrix{Scalar}}
    Qi::Vector{Matrix{Scalar}}

    function SSS{Scalar}(N, n, Gpi, Hpi, Di, Ui, Wi, Vi, Pi, Ri, Qi) where {Scalar<:Number}

        if N < 3
            error("N<3 not supported")
        end

        # check input dimensions
        if length(n) != N || length(Di) != N || length(Ui) != N - 1 ||
           length(Wi) != N - 2 || length(Vi) != N - 1 || length(Pi) != N - 1 ||
           length(Ri) != N - 2 || length(Qi) != N - 1 || length(Hpi) != N - 1 ||
           length(Gpi) != N - 1
            error("Input dimensions do not agree")
        end

        # check diagonal matrix dimensions
        for i = 1:N
            if size(Di[i], 1) != n[i] || size(Di[i], 2) != n[i]
                error("Diagonal matrices dimensions mismatch")
            end
        end

        # check upper triangular matrix dimensions
        for i = 1:N-1
            if size(Ui[i], 1) != n[i] || size(Ui[i], 2) != Hpi[i]
                error("Ui matrices dimensions mismatch")
            end
            if size(Vi[i], 1) != n[i+1] || size(Vi[i], 2) != Hpi[i]
                error("Vi matrices dimensions mismatch")
            end
        end
        for i = 1:N-2
            if size(Wi[i], 1) != Hpi[i] || size(Wi[i], 2) != Hpi[i+1]
                error("Wi matrices dimensions mismatch")
            end
        end

        # check lower triangular matrix dimensions
        for i = 1:N-1
            if size(Pi[i], 1) != n[i+1] || size(Pi[i], 2) != Gpi[i]
                error("Pi matrices dimensions mismatch")
            end
            if size(Qi[i], 1) != n[i] || size(Qi[i], 2) != Gpi[i]
                error("Qi matrices dimensions mismatch")
            end
        end
        for i = 1:N-2
            if size(Ri[i], 1) != Gpi[i+1] || size(Ri[i], 2) != Gpi[i]
                error("Ri translation operators dimensions mismatch")
            end
        end

        new{Scalar}(N, n, Gpi, Hpi,
            map(x -> convert(Matrix{Scalar}, x), Di),
            map(x -> convert(Matrix{Scalar}, x), Ui),
            map(x -> convert(Matrix{Scalar}, x), Wi),
            map(x -> convert(Matrix{Scalar}, x), Vi),
            map(x -> convert(Matrix{Scalar}, x), Pi),
            map(x -> convert(Matrix{Scalar}, x), Ri),
            map(x -> convert(Matrix{Scalar}, x), Qi))


    end

end
Base.:size(A::SSS) = (sum(A.n), sum(A.n))



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











