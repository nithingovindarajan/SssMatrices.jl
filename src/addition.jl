
function add_diagonals_parts(diag1, diag2)
    return [D1i + D2i for (D1i, D2i) ∈ zip(diag1.D, diag2.D)]
end

function add_diagonals_parts_adjoint(diag1, diag2)
    return [D1i + D2i' for (D1i, D2i) ∈ zip(diag1.D, diag2.D)]
end

function add_triangular_parts(triang1, triang2)
    inp = [hcat(inp1i, inp2i) for (inp1i, inp2i) ∈ zip(triang1.inp, triang2.inp)]
    trans = [
        BlockDiagonal([trans1i, trans2i]) for
        (trans1i, trans2i) ∈ zip(triang1.trans, triang2.trans)
    ]
    out = [hcat(out1i, out2i) for (out1i, out2i) ∈ zip(triang1.out, triang2.out)]
    return inp, trans, out
end

function Base.:+(A::SSS, B::SSS)
    @assert A.N == B.N "A and B are not of the same size"
    @assert A.n == B.n "A and B are not conformly partitioned"
    # add diagonal parts
    D = add_diagonals_parts(A.diagonal, B.diagonal)
    Q, R, P = add_triangular_parts(A.lower, B.lower)
    U, W, V = add_triangular_parts(A.upper, B.upper)
    return SSS{Float64}(D, Q, R, P, U, W, V)
end
