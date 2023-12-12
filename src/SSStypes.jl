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
        out::Vector{<:AbstractMatrix{T}},
    ) where {T<:Number}
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
    tot_ranks_upper::Int
    tot_ranks_lower::Int
    off::Vector{Int}
    off_upper::Vector{Int}
    off_lower::Vector{Int}
    block::Vector{UnitRange}
    rblock_upper::Vector{UnitRange}
    rblock_lower::Vector{UnitRange}

    function SSS{T}(D, Q, R, P, U, W, V) where {T<:Number}
        if any([
            !(size(Di, 1) == size(Qi, 1) == size(Pi, 1) == size(Ui, 1) == size(Vi, 1)) for
            (Di, Pi, Qi, Ui, Vi) in zip(D, P, Q, U, V)
        ])
            error("dimensions inconsistent.")
        end

        diagonal = DiagonalPart(convert(Vector{Matrix{T}}, D))
        lower = TriangularPart(
            convert(Vector{Matrix{T}}, Q),
            convert(Vector{Matrix{T}}, R),
            convert(Vector{Matrix{T}}, P),
        )
        upper = TriangularPart(
            convert(Vector{Matrix{T}}, U),
            convert(Vector{Matrix{T}}, W),
            convert(Vector{Matrix{T}}, V),
        )

        n = [size(Di, 1) for Di ∈ D]
        no_blocks = length(D)
        N = sum(n)
        lower_hankel_ranks = [size(Q[i], 2) for i = 1:no_blocks]
        upper_hankel_ranks = [size(V[i], 2) for i = 1:no_blocks]
        tot_ranks_upper = sum(upper_hankel_ranks)
        tot_ranks_lower = sum(lower_hankel_ranks)
        off = [0; cumsum(n)]
        off_upper = [0; cumsum(upper_hankel_ranks)]
        off_lower = [0; cumsum(lower_hankel_ranks)]
        block = [off[i]+1:off[i+1] for i = 1:no_blocks]
        rblock_upper = [off_upper[i]+1:off_upper[i+1] for i = 1:no_blocks]
        rblock_lower = [off_lower[i]+1:off_lower[i+1] for i = 1:no_blocks]

        new{T}(
            diagonal,
            lower,
            upper,
            n,
            no_blocks,
            N,
            lower_hankel_ranks,
            upper_hankel_ranks,
            tot_ranks_upper,
            tot_ranks_lower,
            off,
            off_upper,
            off_lower,
            block,
            rblock_upper,
            rblock_lower,
        )
    end
end
Base.:size(A::SSS) = (A.N, A.N)

function random_SSS(n, ranks_l, ranks_u, seed = 1234)
    @assert length(n) - 1 == length(ranks_l) == length(ranks_u) "input dimensions not compatible"
    Random.seed!(seed)
    no_blocks = length(n)
    D = [rand(n[i], n[i]) for i = 1:no_blocks]
    Q = [
        [rand(n[i], ranks_l[i]) for i = 1:no_blocks-1]
        [rand(n[no_blocks], 0)]
    ]
    R = [
        [rand(ranks_l[1], 0)]
        [rand(ranks_l[i+1], ranks_l[i]) for i = 1:no_blocks-2]
        [rand(0, ranks_l[no_blocks-1])]
    ]
    P = [
        [rand(n[1], 0)]
        [rand(n[i], ranks_l[i-1]) for i = 2:no_blocks]
    ]
    U = [
        [rand(n[i], ranks_u[i]) for i = 1:no_blocks-1]
        [rand(n[no_blocks], 0)]
    ]
    W = [
        [rand(ranks_u[1], 0)]
        [rand(ranks_u[i+1], ranks_u[i]) for i = 1:no_blocks-2]
        [rand(0, ranks_u[no_blocks-1])]
    ]
    V = [
        [rand(n[1], 0)]
        [rand(n[i], ranks_u[i-1]) for i = 2:no_blocks]
    ]
    return SSS{Float64}(D, Q, R, P, U, W, V)
end
