
function add_diagonals_parts(diag1::DiagonalPart, diag2::DiagonalPart)

    return [D1i + D2i for (D1i, D2i) ∈ zip(diag1.D, diag2.D)]


end

function add_diagonals_parts_adjoint(diag1::DiagonalPart, diag2::DiagonalPart)

    return [D1i + D2i' for (D1i, D2i) ∈ zip(diag1.D, diag2.D)]

end


function add_triangular_parts(triang1::TriangularPart, triang2::TriangularPart)

    inp = [hcat(inp1i, inp2i) for (inp1i, inp2i) ∈ zip(triang1.inp, triang2.inp)]
    trans = [BlockDiagonal([trans1i, trans2i]) for (trans1i, trans2i) ∈ zip(triang1.trans, triang2.trans)]
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

    return SSS(D, Q, R, P, U, W, V)

end

function Base.:+(A::SSS, B::Adjoint{Scalar,SSS{Scalar}}) where {Scalar<:Number}

    @assert A.N == B.parent.N "A and B are not of the same size"
    @assert A.n == B.parent.n "A and B are not conformly partitioned"

    # add diagonal parts
    D = add_diagonals_parts_adjoint(A.diagonal, B.parent.diagonal)
    Q, R, P = add_triangular_parts(A.lower, B.parent.upper)
    U, W, V = add_triangular_parts(A.upper, B.parent.lower)

    return SSS(D, Q, R, P, U, W, V)

end
Base.:+(A::Adjoint{Scalar,SSS{Scalar}}, B::SSS) where {Scalar<:Number} = Base.:+(B, A)

function Base.:+(A::Adjoint{ScalarA,SSS{ScalarA}}, B::Adjoint{ScalarB,SSS{ScalarB}}) where {ScalarA<:Number,ScalarB<:Number}
    return (A' + B')'
end