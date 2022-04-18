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

function Sblock(A, Psi, Phi, k)

    return [Matrix(I, size(A.upper.out[k], 2), size(A.upper.out[k], 2)) zeros(size(A.upper.out[k], 2), size(A.lower.inp[k], 2)) -A.upper.out[k]'
        Psi[k] Matrix(I, size(A.lower.inp[k], 2), size(A.lower.inp[k], 2)) -A.lower.inp[k]'
        Phi[k] zeros(A.n[k], size(A.lower.inp[k], 2)) A.diagonal.D[k]]

end

function compute_Phi_Psi(A, Psi, Phi, k)

    a1, a2 = compute_Phi_Psi2(A.upper.out[k-1], A.lower.inp[k-1], A.diagonal.D[k-1], Psi[k-1], Phi[k-1], A.lower.trans[k], A.lower.out[k], A.upper.trans[k-1]', A.upper.inp[k-1])

    return a1, a2

end


function compute_Phi_Psi2(V_i,
    Q_i,
    D_i,
    Psi_i,
    Phi_i,
    R_iplus1,
    P_iplus1,
    W_i,
    U_i)


    Temp = (Psi_i * W_i - (-Q_i' + Psi_i * V_i') * inv(D_i + Phi_i * V_i') * (Phi_i * W_i + U_i))

    return R_iplus1 * Temp, -P_iplus1 * Temp

end


function schurcomplement(A11, A12, A21, A22)

    return A22 - A21 * (A11 \ A12)

end


function compute_schurcomplements(A)

    Psi = SchurComplements(undef, A.no_blocks)
    Phi = SchurComplements(undef, A.no_blocks)

    Psi[1] = zeros(size(A.lower.inp[1], 2), size(A.upper.out[1], 2))
    Phi[1] = zeros(A.n[1], size(A.upper.out[1], 2))

    for k = 2:A.no_blocks

        Psi[k], Phi[k] = compute_Phi_Psi(A, Psi, Phi, k)

    end

    return Psi, Phi

end

function forwardsolve(A, k, Psi, Phi, c, q, r, s)

    c_lift = [zeros(A.upper_hankel_ranks[k])
        zeros(A.lower_hankel_ranks[k])
        c[A.block[k]]]



    x_lift = c_lift - subdiagblock(A, k) * (Sblock(A, Psi, Phi, k - 1) \ [q[A.rblock_upper[k-1]]
        r[A.rblock_lower[k-1]]
        s[A.block[k-1]]])

    return x_lift[1:A.upper_hankel_ranks[k]], x_lift[(A.upper_hankel_ranks[k]+1):A.upper_hankel_ranks[k]+A.lower_hankel_ranks[k]], x_lift[(A.upper_hankel_ranks[k]+A.lower_hankel_ranks[k]+1):end]
end

function forwarditerate(A, Psi, Phi, c)


    q = zeros(ComplexF64, A.tot_ranks_upper)
    r = zeros(ComplexF64, A.tot_ranks_lower)
    s = zeros(ComplexF64, A.N)

    s[A.block[1]] = c[A.block[1]]

    for k = 2:A.no_blocks

        q[A.rblock_upper[k]], r[A.rblock_lower[k]], s[A.block[k]] = forwardsolve(A, k, Psi, Phi, c, q, r, s)

    end


    return q, r, s

end



function backwarditerate!(q, r, s, A, Psi, Phi)




    x_lift = Sblock(A, Psi, Phi, A.no_blocks) \ [q[A.rblock_upper[A.no_blocks]]
        r[A.rblock_lower[A.no_blocks]]
        s[A.block[A.no_blocks]]]

    q[A.rblock_upper[A.no_blocks]] = x_lift[1:A.upper_hankel_ranks[A.no_blocks]]
    r[A.rblock_lower[A.no_blocks]] = x_lift[A.upper_hankel_ranks[A.no_blocks]+1:A.lower_hankel_ranks[A.no_blocks]+A.upper_hankel_ranks[A.no_blocks]]
    s[A.block[A.no_blocks]] = x_lift[A.lower_hankel_ranks[A.no_blocks]+A.upper_hankel_ranks[A.no_blocks]+1:end]

    for k = A.no_blocks:-1:2


        x_lift = Sblock(A, Psi, Phi, k - 1) \ ([q[A.rblock_upper[k-1]]
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
        Psi, Phi = compute_schurcomplements(A)
        q, r, s = forwarditerate(A, Psi, Phi, c)
        backwarditerate!(q, r, s, A, Psi, Phi)

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


n = 5
r = 3
V_i = rand(n, r)
Q_i = rand(n, r)
D_i = rand(n, n)
Psi_i = rand(n, r)
Phi_i = rand(n, r)
R_iplus1 = rand(r, r)
P_iplus1 = rand(n, r)
W_i = rand(r, r)
U_i = rand(n, r)


V_iplus1 = rand(n, r)
Q_iplus1 = rand(n, r)
D_iplus1 = rand(n, n)


Psi_i = rand(r, r)
Phi_i = rand(n, r)




A = [Matrix(I, r, r) zeros(r, r) -V_iplus1'
    zeros(r, r) Matrix(I, r, r) -Q_iplus1'
    zeros(n, r) zeros(n, r) D_iplus1]
B = [zeros(r, r) zeros(r, r) zeros(r, n)
    zeros(r, r) -R_iplus1 zeros(r, n)
    zeros(n, r) P_iplus1 zeros(n, n)]
D = [Matrix(I, r, r) zeros(r, r) -V_i'
    Psi_i Matrix(I, r, r) -Q_i'
    Phi_i zeros(n, r) D_i]
C = [-W_i zeros(r, r) zeros(r, n)
    zeros(r, r) zeros(r, r) zeros(r, n)
    U_i zeros(n, r) zeros(n, n)]

S = A - B * inv(D) * C


# iteration
Psi_iplus1, Phi_iplus1 = compute_Phi_Psi2(V_i, Q_i, D_i, Psi_i, Phi_i, R_iplus1, P_iplus1, W_i, U_i)

