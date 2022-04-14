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


function forwarditerate(A, S, c)

    x_lift = Vector{Vector}(undef, A.no_blocks)

    x_lift[1] = [zeros(A.upper_hankel_ranks[1]); zeros(A.lower_hankel_ranks[1]); c[A.off[1]+1:A.off[2]]]

    for k = 2:A.no_blocks

        c_lift = [zeros(A.upper_hankel_ranks[k])
            zeros(A.lower_hankel_ranks[k])
            c[A.off[k]+1:A.off[k+1]]]


        x_lift[k] = c_lift - subdiagblock(A, k) * (S[k-1] \ x_lift[k-1])

    end


    return x_lift

end



function backwarditerate!(x_lift, A, S)


    x_lift[A.no_blocks] = S[A.no_blocks] \ x_lift[A.no_blocks]

    for k = A.no_blocks:-1:2

        x_lift[k-1] = S[k-1] \ (x_lift[k-1] - supdiagblock(A, k) * x_lift[k])

    end


end



function Base.:\(A::SSS, c::Vector)



    if A.no_blocks == 1
        return A.diagonal.D[1] \ c
    else

        # note: D^-1 is computed twice. It might be good to save the QR / LU of Di in memory?
        S = compute_schurcomplements(A)
        x_lift = forwarditerate(A, S, c)
        t = deepcopy(x_lift)
        backwarditerate!(x_lift, A, S)

        x = reduce(vcat, (x_lift[i][A.upper_hankel_ranks[i]+A.lower_hankel_ranks[i]+1:end] for i = 1:A.no_blocks))

        return x

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














# Experiment

# generators
n = 4
r = 2
D1, D2, D3 = rand(n, n), rand(n, n), rand(n, n)
Q1, Q2, Q3 = rand(n, r), rand(n, 2), rand(n, 0)
R1, R2, R3 = rand(r, 0), rand(r, r), rand(0, r)
P1, P2, P3 = rand(n, 0), rand(n, r), rand(n, r)
U1, U2, U3 = rand(n, r), rand(n, r), rand(n, 0)
W1, W2, W3 = rand(0, r), rand(r, r), rand(r, 0)
V1, V2, V3 = rand(n, 0), rand(n, r), rand(n, r)

# right hand side
c1, c2, c3 = rand(n), rand(n), rand(n)

# The SSS matrix
A_SSS = [D1 U1*V2' U1*W2*V3'
    P2*Q1' D2 U2*V3'
    P3*R2*Q1' P3*Q2' D3]

A11 = [Matrix(I, 0, 0) zeros(0, r) -V1'
    zeros(r, 0) Matrix(I, r, r) -Q1'
    zeros(n, 0) zeros(n, r) D1]


A21 = [zeros(r, 0) zeros(r, r) zeros(r, n)
    zeros(r, 0) -R2 zeros(r, n)
    zeros(n, 0) P2 zeros(n, n)]

A31 = [zeros(r, 0) zeros(r, r) zeros(r, n)
    zeros(0, 0) zeros(0, r) zeros(0, n)
    zeros(n, 0) zeros(n, r) zeros(n, n)]

A12 = [-W1 zeros(0, 2) zeros(0, n)
    zeros(r, r) zeros(r, 2) zeros(2, n)
    U1 zeros(n, 2) zeros(n, n)]

A22 = [Matrix(I, r, r) zeros(r, r) -V2'
    zeros(r, r) Matrix(I, r, r) -Q2'
    zeros(n, r) zeros(n, r) D2]

A32 = [zeros(r, r) zeros(r, r) zeros(r, n)
    zeros(0, r) -R3 zeros(0, n)
    zeros(n, r) P3 zeros(n, n)]

A13 = [zeros(0, r) zeros(0, 0) zeros(0, n)
    zeros(r, r) zeros(r, 0) zeros(r, n)
    zeros(n, r) zeros(n, 0) zeros(n, n)]

A23 = [-W2 zeros(r, 0) zeros(r, n)
    zeros(r, r) zeros(r, 0) zeros(r, n)
    U2 zeros(n, 0) zeros(n, n)]

A33 = [Matrix(I, r, r) zeros(r, 0) -V3'
    zeros(0, r) Matrix(I, 0, 0) -Q3'
    zeros(n, r) zeros(n, 0) D3]


Asparse = [A11 A12 A13
    A21 A22 A23
    A31 A32 A33]


#test if embedding is constructed properly
x_ref = A_SSS \ [c1; c2; c3]
largevec = Asparse \ [zeros(0); zeros(r); c1; zeros(r); zeros(r); c2; zeros(r); zeros(0); c3]
x = [largevec[(0+r+1):(0+r+n)]; largevec[(0+r+n+r+r+1):(0+r+n+r+r+n)]; largevec[(0+r+n+r+n+r+r+0+1):end]]
x ≈ x_ref

