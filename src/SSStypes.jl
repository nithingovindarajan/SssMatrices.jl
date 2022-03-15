export SSSFamily, DiagonalSSS, StrictlyLowerTriangularSSS, StrictlyUpperTriangularSSS, LowerTriangularSSS, UpperTriangularSSS, SSS


abstract type SSSFamily{Scalar<:Number} <: AbstractMatrix{Scalar} end


struct DiagonalSSS{Scalar<:Number} <: SSSFamily{Scalar}

    # attributes required for construction
    D::Vector{AbstractMatrix{Scalar}}

    # attributes generated in construction
    n::Vector{Int64}
    no_blocks::Int64
    N::Int64


    function DiagonalSSS{Scalar}(D) where {Scalar<:Number}

        if any([size(Di, 1) != size(Di, 2) for Di in D])
            error("All diagonal block entries must be square.")
        end

        n = [size(Di, 1) for Di in D]
        no_blocks = length(D)
        N = sum(n)


        new{Scalar}(map(x -> convert(Matrix{Scalar}, x), D), n, no_blocks, N)
    end

end
Base.:size(A::DiagonalSSS) = (A.N, A.N)




struct TriangularSSS{Scalar<:Number}

    # attributes required for construction
    inp::Vector{AbstractMatrix{Scalar}}
    trans::Vector{AbstractMatrix{Scalar}}
    out::Vector{AbstractMatrix{Scalar}}

    # attributes generated in construction
    n::Vector{Int64}
    no_blocks::Int64
    N::Int64
    hankel_ranks::Vector{Int64}   # hankel block ranks lower triangular part


    function TriangularSSS{Scalar}(inp, trans, out) where {Scalar<:Number}


        # check input dimensions
        if !(length(out) == length(inp) == length(trans))
            error("Input dimensions do not agree")
        end
        no_blocks = length(out)

        if any([size(out_i, 1) != size(inp_i, 1) for (out_i, inp_i) in zip(out, inp)])
            error("dimensions inconsistent.")
        end
        n = [size(out_i, 1) for out_i in out]
        N = sum(n)

        # check lower triangular matrix dimensions
        for i = 1:no_blocks-1
            if size(out[i+1], 2) != size(inp[i], 2)
                error("P and Q matrices dimensions mismatch")
            end
        end
        for i = 1:no_blocks-2
            if size(out[i], 1) != size(trans[i+1], 2)
                error("Q and R matrices dimensions mismatch")
            end
        end
        for i = 2:no_blocks-2
            if size(trans[i], 1) != size(trans[i+1], 2)
                error("R translation operators dimensions mismatch")
            end
        end
        for i = 3:no_blocks
            if size(inp[i], 2) != size(trans[i-1], 1)
                error("R translation operators dimensions mismatch")
            end
        end

        hankel_ranks = [size(inp[i], 2) for i in 1:(no_blocks-1)]

        new{Scalar}(map(x -> convert(Matrix{Scalar}, x), inp),
            map(x -> convert(Matrix{Scalar}, x), trans),
            map(x -> convert(Matrix{Scalar}, x), out),
            N, n, no_blocks, hankel_ranks)

    end

end




struct StrictlyLowerTriangularSSS{Scalar<:Number} <: SSSFamily{Scalar}


    # attributes generated in construction
    lower::TriangularSSS{Scalar}

    function StrictlyLowerTriangularSSS{Scalar}(Q, R, P) where {Scalar<:Number}

        new{Scalar}(TriangularSSS{Scalar}(Q, R, P))

    end

end
Base.:size(A::StrictlyLowerTriangularSSS) = (A.lower.N, A.lower.N)




struct StrictlyUpperTriangularSSS{Scalar<:Number} <: AbstractMatrix{Scalar}

    # attributes generated in construction
    upper::TriangularSSS{Scalar}

    function StrictlyUpperTriangularSSS{Scalar}(V, W, U) where {Scalar<:Number}

        new{Scalar}(TriangularSSS{Scalar}(V, W, U))

    end

end
Base.:size(A::StrictlyUpperTriangularSSS) = (A.upper.N, A.upper.N)




struct UpperTriangularSSS{Scalar<:Number} <: AbstractMatrix{Scalar}

    # attributes generated in construction
    diagonal::DiagonalSSS{Scalar}
    upper::StrictlyUpperTriangularSSS{Scalar}

    function UpperTriangularSSS{Scalar}(D, Q, R, P, V, W, U) where {Scalar<:Number}

        diagonal = DiagonalSSS{Scalar}(D)
        upper = StrictlyUpperTriangularSSS{Scalar}(V, W, U)


        if !(diagonal.no_blocks == upper.no_blocks)
            error("No blocks not consistent")
        end
        if any(diagonal.n .!= upper.n)
            error("dimensions of blocks not consistent")
        end


        new{Scalar}(diagonal, upper)

    end

end
Base.:size(A::UpperTriangularSSS) = (A.diagonal.N, A.diagonal.N)


struct LowerTriangularSSS{Scalar<:Number} <: AbstractMatrix{Scalar}

    # attributes generated in construction
    diagonal::DiagonalSSS{Scalar}
    lower::StrictlyLowerTriangularSSS{Scalar}

    function LowerTriangularSSS{Scalar}(D, Q, R, P, V, W, U) where {Scalar<:Number}

        diagonal = DiagonalSSS{Scalar}(D)
        lower = StrictlyLowerTriangularSSS{Scalar}(V, W, U)


        if !(diagonal.no_blocks == lower.no_blocks)
            error("No blocks not consistent")
        end
        if any(diagonal.n .!= lower.n)
            error("dimensions of blocks not consistent")
        end


        new{Scalar}(diagonal, lower)

    end

end
Base.:size(A::LowerTriangularSSS) = (A.diagonal.N, A.diagonal.N)




struct SSS{Scalar<:Number} <: AbstractMatrix{Scalar}

    # attributes generated in construction
    diagonal::DiagonalSSS{Scalar}
    lower::StrictlyLowerTriangularSSS{Scalar}
    upper::StrictlyUpperTriangularSSS{Scalar}

    function SSS{Scalar}(D, Q, R, P, V, W, U) where {Scalar<:Number}

        diagonal = DiagonalSSS{Scalar}(D)
        lower = StrictlyLowerTriangularSSS{Scalar}(R, P, V)
        upper = StrictlyUpperTriangularSSS{Scalar}(V, W, U)


        if !(diagonal.no_blocks == lower.no_blocks == upper.no_blocks)
            error("No blocks not consistent")
        end
        if any(diagonal.n .!= lower.n .!= upper.n)
            error("dimensions of blocks not consistent")
        end


        new{Scalar}(diagonal, lower, upper)

    end

end
Base.:size(A::SSS) = (A.diagonal.N, A.diagonal.N)

