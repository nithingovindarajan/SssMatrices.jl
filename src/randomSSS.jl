########### random SSS matrix generation ############

export random_SSS

function random_SSS(n::Vector{<:Integer}, ranks_l::Vector{<:Integer}, ranks_u::Vector{<:Integer})

    @assert length(n) - 1 == length(ranks_l) == length(ranks_u) "input dimensions not compatible"


    no_blocks = length(n)
    D = [rand(n[i], n[i]) for i = 1:no_blocks]

    Q = [[rand(n[i], ranks_l[i]) for i = 1:no_blocks-1]
        [rand(n[no_blocks], 0)]]
    R = [[rand(ranks_l[1], 0)]
        [rand(ranks_l[i+1], ranks_l[i]) for i = 1:no_blocks-2]
        [rand(0, ranks_l[no_blocks-1])]]
    P = [[rand(n[1], 0)]
        [rand(n[i], ranks_l[i-1]) for i = 2:no_blocks]]

    U = [[rand(n[i], ranks_u[i]) for i = 1:no_blocks-1]
        [rand(n[no_blocks], 0)]]
    W = [[rand(ranks_u[1], 0)]
        [rand(ranks_u[i+1], ranks_u[i]) for i = 1:no_blocks-2]
        [rand(0, ranks_u[no_blocks-1])]]
    V = [[rand(n[1], 0)]
        [rand(n[i], ranks_u[i-1]) for i = 2:no_blocks]]

    return SSS(D, Q, R, P, U, W, V)

end
