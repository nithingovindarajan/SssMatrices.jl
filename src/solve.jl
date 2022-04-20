
## All them of are computing inverse of (D_i + Phi_i * V_i') !!
function compute_Phi_Psi(V_i,
    Q_i,
    D_i,
    Psi_i,
    Phi_i,
    R_iplus1,
    P_iplus1,
    W_i,
    U_i)


    temp = Psi_i * W_i - (-Q_i' + Psi_i * V_i') * ((D_i + Phi_i * V_i') \ (Phi_i * W_i + U_i))

    return R_iplus1 * temp, -P_iplus1 * temp

end

function forwardsolve(D_kmin1, Q_kmin1, V_kmin1, R_k, P_k, Psi_kmin1, Phi_kmin1, c_k, r_kmin1, s_kmin1)


    temp = r_kmin1 - (-Q_kmin1' + Psi_kmin1 * V_kmin1') * ((D_kmin1 + Phi_kmin1 * V_kmin1') \ s_kmin1)

    return R_k * temp, c_k - P_k * temp
end


function backwardsolve(D_kmin1, V_kmin1, Phi_kmin1, q, s)


    xk = (D_kmin1 + Phi_kmin1 * V_kmin1') \ (-Phi_kmin1 * q + s)


    return xk, q + V_kmin1' * xk
end


export compute_schurcomplements

function compute_schurcomplements(A)

    Psi = Vector{Matrix}(undef, A.no_blocks)
    Phi = Vector{Matrix}(undef, A.no_blocks)

    Psi[1] = zeros(size(A.lower.inp[1], 2), size(A.upper.out[1], 2))
    Phi[1] = zeros(A.n[1], size(A.upper.out[1], 2))

    for k = 2:A.no_blocks

        Psi[k], Phi[k] = compute_Phi_Psi(A.upper.out[k-1], A.lower.inp[k-1], A.diagonal.D[k-1], Psi[k-1], Phi[k-1],
            A.lower.trans[k], A.lower.out[k], A.upper.trans[k-1]', A.upper.inp[k-1])

    end

    return Psi, Phi

end



function forwarditerate(A, Psi, Phi, c)


    r = zeros(A.lower_hankel_ranks[1])
    x = zeros(ComplexF64, A.N)

    x[A.block[1]] = c[A.block[1]]

    for k = 2:A.no_blocks

        r, x[A.block[k]] = forwardsolve(A.diagonal.D[k-1], A.lower.inp[k-1], A.upper.out[k-1],
            A.lower.trans[k], A.lower.out[k], Psi[k-1], Phi[k-1], c[A.block[k]], r, x[A.block[k-1]])
    end


    return x

end



function backwarditerate!(x, A, Phi)


    x[A.block[A.no_blocks]], q = backwardsolve(A.diagonal.D[A.no_blocks], A.upper.out[A.no_blocks], Phi[A.no_blocks], zeros(A.upper_hankel_ranks[A.no_blocks]), x[A.block[A.no_blocks]])

    for k = A.no_blocks:-1:2

        x[A.block[k-1]], q = backwardsolve(A.diagonal.D[k-1], A.upper.out[k-1], Phi[k-1], A.upper.trans[k-1]' * q, x[A.block[k-1]] - A.upper.inp[k-1] * q)
    end

end



function Base.:\(A::SSS, c::Vector)



    if A.no_blocks == 1
        return A.diagonal.D[1] \ c
    else

        # note: D^-1 is computed twice. It might be good to save the QR / LU of Di in memory?
        Psi, Phi = compute_schurcomplements(A)
        x = forwarditerate(A, Psi, Phi, c)
        backwarditerate!(x, A, Phi)

        return x

    end

end



