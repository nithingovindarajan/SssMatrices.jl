function lowrankapprox(B::AbstractMatrix, threshold::Float64)

    # compute SVD
    U, sigma, V = svd(B)

    # rank
    p = findlast(x -> x > sigma[1] * threshold, sigma)
    if p == nothing
        p = 0
    end

    #truncate
    X = U[:, 1:p]
    Y = V[:, 1:p] * Diagonal(sigma[1:p])

    return X, Y
end



### Iterating through blocks ###

# struct Blocks
#     n::Vector{Int}
#     start_pt::Int
#     end_pt::Int


#     function Blocks(n::Vector{Int}, start_pt::Int, end_pt::Int)

#         @assert 1 <= start_pt <= end_pt <= length(n)

#         new(n, start_pt, end_pt)
#     end
# end
# Blocks(n::Vector{Int}) = Blocks(n, 1, length(n))
# Blocks(n::Vector{Int}, start_pt::Int, end_pt::Nothing) = Blocks(n, start_pt, length(n))
# Blocks(n::Vector{Int}, start_pt::Nothing, end_pt::Int) = Blocks(n, 1, end_pt)



# function Base.iterate(blocks::Blocks)

#     if isempty(blocks.n)
#         return nothing

#     else

#         entry = blocks.start_pt
#         off = blocks.n[entry]
#         item = 1:off

#         return (entry, item), (entry, off)
#     end


# end


# function Base.iterate(blocks::BlockSizes, state)

#     if state[1] == blocks.end_pt

#         return nothing

#     else

#         entry = state[1] + 1
#         off = state[2] + blocks.n[entry]
#         item = (state[2]+1):off



#         return (entry, item), (entry, off)

#     end



# end


