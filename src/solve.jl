export forwarditerate

SchurComplements = Vector{Matrix}

function diagblock(A, k)

    return [Matrix(I, size(A.upper.out[k], 2), size(A.upper.out[k], 2)) zeros(size(A.upper.out[k], 2), size(A.lower.inp[k], 2)) -A.upper.out[k]'
        zeros(size(A.lower.inp[k], 2), size(A.upper.out[k], 2)) Matrix(I, size(A.lower.inp[k], 2), size(A.lower.inp[k], 2)) -A.lower.inp[k]'
        zeros(A.n[k], size(A.upper.out[k], 2)) zeros(A.n[k], size(A.lower.inp[k], 2)) A.diagonal.D[k]]

end

function supdiagblock(A, k)
    return [-A.upper.trans[k-1]' zeros(size(A.upper.out[k-1], 2), size(A.lower.inp[k], 2)) zeros(size(A.upper.out[k-1], 2), A.n[k])
        zeros(size(A.lower.inp[k-1], 2), size(A.upper.inp[k-1], 2)) zeros(size(A.lower.inp[k-1], 2), size(A.lower.inp[k], 2)) zeros(size(A.lower.inp[k-1], 2), A.n[k])
        A.upper.inp[k-1] zeros(A.n[k-1], size(A.lower.inp[k], 2)) zeros(A.n[k-1], A.n[k])]
end

function subdiagblock(A, k)
    return [zeros(size(A.upper.out[k], 2), size(A.upper.out[k-1], 2)) zeros(size(A.upper.out[k], 2), size(A.lower.out[k], 2)) zeros(size(A.upper.out[k], 2), A.n[k-1])
        zeros(size(A.lower.inp[k], 2), size(A.upper.out[k-1], 2)) -A.lower.trans[k] zeros(size(A.lower.inp[k], 2), A.n[k-1])
        zeros(A.n[k], size(A.upper.out[k-1], 2)) A.lower.out[k] zeros(A.n[k], A.n[k-1])]

end


function schurcomplement(A11, A12, A21, A22)

    return A22 - A21 * (A11 \ A12)

end


function compute_schurcomplements(A)

    S = SchurComplements(undef, A.no_blocks)

    S[1] = diagblock(A, 1)

    for k = 2:A.no_blocks

        S[k] = schurcomplement(S[k-1], supdiagblock(A, k), subdiagblock(A, k), diagblock(A, k))

    end

    return S

end

function forwardsolve(A, k, S, c, q, r, s)

    c_lift = [zeros(A.upper_hankel_ranks[k])
        zeros(A.lower_hankel_ranks[k])
        c[A.block[k]]]

    x_lift = c_lift - subdiagblock(A, k) * (S[k-1] \ [q[A.rblock_upper[k-1]]
        r[A.rblock_lower[k-1]]
        s[A.block[k-1]]])

    return x_lift[1:A.upper_hankel_ranks[k]], x_lift[(A.upper_hankel_ranks[k]+1):A.upper_hankel_ranks[k]+A.lower_hankel_ranks[k]], x_lift[(A.upper_hankel_ranks[k]+A.lower_hankel_ranks[k]+1):end]
end

function forwarditerate(A, S, c)


    q = zeros(ComplexF64, A.tot_ranks_upper)
    r = zeros(ComplexF64, A.tot_ranks_lower)
    s = zeros(ComplexF64, A.N)

    s[A.block[1]] = c[A.block[1]]

    for k = 2:A.no_blocks

        q[A.rblock_upper[k]], r[A.rblock_lower[k]], s[A.block[k]] = forwardsolve(A, k, S, c, q, r, s)

    end


    return q, r, s

end



function backwarditerate!(q, r, s, A, S)




    x_lift = S[A.no_blocks] \ [q[A.rblock_upper[A.no_blocks]]
        r[A.rblock_lower[A.no_blocks]]
        s[A.block[A.no_blocks]]]

    q[A.rblock_upper[A.no_blocks]] = x_lift[1:A.upper_hankel_ranks[A.no_blocks]]
    r[A.rblock_lower[A.no_blocks]] = x_lift[A.upper_hankel_ranks[A.no_blocks]+1:A.lower_hankel_ranks[A.no_blocks]+A.upper_hankel_ranks[A.no_blocks]]
    s[A.block[A.no_blocks]] = x_lift[A.lower_hankel_ranks[A.no_blocks]+A.upper_hankel_ranks[A.no_blocks]+1:end]

    for k = A.no_blocks:-1:2


        x_lift = S[k-1] \ ([q[A.rblock_upper[k-1]]
            r[A.rblock_lower[k-1]]
            s[A.block[k-1]]] - supdiagblock(A, k) * [q[A.rblock_upper[k]]
            r[A.rblock_lower[k]]
            s[A.block[k]]])

        q[A.rblock_upper[k-1]] = x_lift[1:A.upper_hankel_ranks[k-1]]
        r[A.rblock_lower[k-1]] = x_lift[A.upper_hankel_ranks[k-1]+1:A.lower_hankel_ranks[k-1]+A.upper_hankel_ranks[k-1]]
        s[A.block[k-1]] = x_lift[A.lower_hankel_ranks[k-1]+A.upper_hankel_ranks[k-1]+1:end]

    end

end



function Base.:\(A::SSS, c::Vector)



    if A.no_blocks == 1
        return A.diagonal.D[1] \ c
    else

        # note: D^-1 is computed twice. It might be good to save the QR / LU of Di in memory?
        S = compute_schurcomplements(A)
        q, r, s = forwarditerate(A, S, c)
        backwarditerate!(q, r, s, A, S)

        return s

    end

end


# function Base.:\(A::SSS, S::SchurComplements, c::Vector)



#     if A.no_blocks == 1
#         return A.diagonal.D[1] \ c
#     else


#         x_lift = forwarditerate(A, S, c)
#         backwarditerate!(x_lift, A, S)

#         return x

#     end

# end








