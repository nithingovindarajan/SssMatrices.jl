struct DiagonalPart{Scalar<:Number}

    # attributes required for construction
    D::Vector{AbstractMatrix{Scalar}}

    function DiagonalPart(D::Vector{<:AbstractMatrix{T}}) where {T<:Number}

        if any([size(Di, 1) != size(Di, 2) for Di ∈ D])
            error("All diagonal block entries must be square.")
        end
        Scalar = eltype(eltype(D))
        new{Scalar}(D)

    end

end

struct TriangularPart{Scalar<:Number}

    # attributes required for construction
    inp::Vector{AbstractMatrix{Scalar}}
    trans::Vector{AbstractMatrix{Scalar}}
    out::Vector{AbstractMatrix{Scalar}}



    function TriangularPart(
        inp::Vector{<:AbstractMatrix{T}},
        trans::Vector{<:AbstractMatrix{T}},
        out::Vector{<:AbstractMatrix{T}}) where {T<:Number}


        # check input dimensions
        if !(length(out) == length(inp) == length(trans))
            error("Input dimensions do not agree")
        end
        no_blocks = length(out)

        if any([size(out_i, 1) != size(inp_i, 1) for (out_i, inp_i) ∈ zip(out, inp)])
            error("dimensions inconsistent.")
        end

        #TODO write this more efficiently...

        # check lower triangular matrix dimensions

        for i = 1:no_blocks-1
            if size(out[i+1], 2) != size(inp[i], 2)
                error("P and Q matrices dimensions mismatch")
            end
        end
        for i = 1:no_blocks-2
            if size(inp[i], 2) != size(trans[i+1], 2)
                error("Q and R matrices dimensions mismatch")
            end
        end
        for i = 1:no_blocks-1
            if size(trans[i], 1) != size(trans[i+1], 2)
                error("R translation operators dimensions mismatch")
            end
        end
        for i = 3:no_blocks
            if size(out[i], 2) != size(trans[i-1], 1)
                error("R translation operators dimensions mismatch")
            end
        end
        @assert size(inp[end], 2) == 0
        @assert size(out[1], 2) == 0
        @assert size(trans[1], 2) == 0
        @assert size(trans[end], 1) == 0


        Scalar = eltype(eltype(inp))
        new{Scalar}(inp, trans, out)

    end

end



########### SSS types ############


export SSSFamily, DiagonalSSS, SSS
# To be done: StrictlyLowerTriangularSSS, StrictlyUpperTriangularSSS, LowerTriangularSSS, UpperTriangularSSS, HermitianSSS

abstract type SSSFamily{Scalar<:Number} <: AbstractMatrix{Scalar} end


struct DiagonalSSS{Scalar<:Number} <: SSSFamily{Scalar}


    # attributes generated in construction
    diagonal::DiagonalPart{Scalar}
    n::Vector{Int}
    no_blocks::Int
    N::Int
    off::Vector{Int}

    function DiagonalSSS(D::Vector{<:AbstractMatrix{T}}) where {T<:Number}

        diagonal = DiagonalPart(D)

        n = [size(Di, 1) for Di ∈ diagonal.D]
        no_blocks = length(diagonal.D)
        N = sum(n)
        off = [0; cumsum(n)]
        new{T}(diagonal, n, no_blocks, N, off)
    end

end



struct SSS{Scalar<:Number} <: AbstractMatrix{Scalar}

    # attributes generated in construction
    diagonal::DiagonalPart{Scalar}
    lower::TriangularPart{Scalar}
    upper::TriangularPart{Scalar}

    n::Vector{Int}
    no_blocks::Int
    N::Int
    lower_hankel_ranks::Vector{Int}
    upper_hankel_ranks::Vector{Int}
    off::Vector{Int}

    function SSS{T}(D, Q, R, P, U, W, V) where {T<:Number}

        if any([!(size(Di, 1) == size(Qi, 1) == size(Pi, 1) == size(Ui, 1) == size(Vi, 1)) for (Di, Pi, Qi, Ui, Vi) in zip(D, P, Q, U, V)])
            error("dimensions inconsistent.")
        end


        diagonal = DiagonalPart(convert(Vector{Matrix{T}}, D))
        lower = TriangularPart(convert(Vector{Matrix{T}}, Q), convert(Vector{Matrix{T}}, R), convert(Vector{Matrix{T}}, P))
        upper = TriangularPart(convert(Vector{Matrix{T}}, U), convert(Vector{Matrix{T}}, W), convert(Vector{Matrix{T}}, V))

        n = [size(Di, 1) for Di ∈ D]
        no_blocks = length(D)
        N = sum(n)
        lower_hankel_ranks = [size(Q[i], 2) for i = 1:no_blocks-1]
        upper_hankel_ranks = [size(U[i], 2) for i = 1:no_blocks-1]
        off = [0; cumsum(n)]

        new{T}(diagonal, lower, upper, n, no_blocks, N, lower_hankel_ranks, upper_hankel_ranks, off)

    end

end




