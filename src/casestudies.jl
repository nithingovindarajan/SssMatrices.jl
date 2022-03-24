export Δ_1d, Δ_2d, Δ_2d_nested_dissection, Δ_1d_schurcompl, semisepbandedmatrix, cauchyish

function Δ_1d(n)

    Is = [1:n; 1:n-1; 2:n]
    Js = [1:n; 2:n; 1:n-1]
    Vs = [fill(2, n); fill(-1, 2n - 2)]

    return sparse(Is, Js, Vs)

end

function Δ_2d(n, m)

    return kron(Δ_1d(n), sparse(I, m, m)) + kron(sparse(I, n, n), Δ_1d(m))

end



function nest(m, n, ndissect=30, lm=m)
    # https://github.com/mohamed82008/18S096-iap17/blob/master/lecture5/Nested-Dissection.ipynb

    if ndissect <= 0 || m * n < 5
        return reshape([i + (j - 1) * lm for i = 1:m, j = 1:n], m * n)
    elseif m >= n
        msep = div(m, 2)
        N1 = nest(msep - 1, n, ndissect - 1, lm)
        N2 = nest(m - msep, n, ndissect - 1, lm) + msep
        Ns = msep + (0:n-1) * lm
        return [N1; N2; Ns]
    else
        nsep = div(n, 2)
        N1 = nest(m, nsep - 1, ndissect - 1, lm)
        N2 = nest(m, n - nsep, ndissect - 1, lm) + nsep * lm
        Ns = (1:m) + (nsep - 1) * lm
        return [N1; N2; Ns]
    end
end


function Δ_2d_nested_dissection(n, m)

    return Δ_2d(n, m)

end


function Δ_1d_schurcompl(n)

    A = Δ_1d(2 * n)
    S = A[1:n, 1:n] - A[1:n, n+1:end] * inv(Matrix(A[n+1:end, n+1:end])) * A[n+1:end, 1:n]

    return S

end




function semisepbandedmatrix(n; r_lower=1, r_upper=1, band_lower=1, band_upper=1, seed=1234)
    # This function creates a

    # set seed
    Random.seed!(seed)

    # intialize
    A = zeros(n, n)

    # semi separable part
    A = A + Array(LowerTriangular(rand(n, r_lower) * rand(r_lower, n))) +
        Array(UpperTriangular(rand(n, r_upper) * rand(r_upper, n)))

    # Banded part
    for l = band_lower:0
        for i = 1:n+l
            A[i-l, i] = A[i-l, i] + rand(Float64)
        end
    end
    for l = 1:band_upper
        for i = 1:n-l
            A[i, i+l] = A[i, i+l] + rand(Float64)
        end
    end

    return A

end

function cauchyish(n)

    θ = collect(range(2 * pi / n, length=n, stop=2 * pi))
    B = zeros(n, n)

    for i = 1:n
        for j = 1:n

            if i != j
                B[i, j] = 1 / abs(exp(im * θ[i]) - exp(im * θ[j]))
            else
                B[i, j] = 0
            end

        end
    end

    return B

end


