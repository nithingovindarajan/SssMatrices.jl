using SparseArrays

function Δ_1d(n)
    Is = [1:n; 1:n-1; 2:n]
    Js = [1:n; 2:n; 1:n-1]
    Vs = [fill(2, n); fill(-1, 2n - 2)]
    return Matrix(sparse(Is, Js, Vs))
end

function Δ_2d(n, m)
    return kron(Δ_1d(n), sparse(I, m, m)) + kron(sparse(I, n, n), Δ_1d(m))
end

function nest_diss_no_sep(off_x, off_y, m, n, M, N; threshold = 2)
    if max(m, n) < threshold
        return [(off_x - 1 + i) + (off_y - 2 + j) * M for j = 1:n for i = 1:m]   # CartesianIndex(off_x - 1 + i, off_y - 1 + j)
    elseif m >= n
        msep = div(m, 2)
        N1 = nest_diss_no_sep(off_x, off_y, msep, n, M, N)
        N2 = nest_diss_no_sep(off_x + msep, off_y, m - msep, n, M, N)
        return [N1; N2]
    else
        nsep = div(n, 2)
        N1 = nest_diss_no_sep(off_x, off_y, m, nsep, M, N)
        N2 = nest_diss_no_sep(off_x, off_y + nsep, m, n - nsep, M, N)
        return [N1; N2]
    end
end
nest_diss_no_sep(M, N; threshold = 2) =
    nest_diss_no_sep(1, 1, M, N, M, N; threshold = threshold)
nest_diss_no_sep(M; threshold = 2) = nest_diss_no_sep(M, M; threshold = threshold)

function Δ_2d_nested_dissection(n, m)
    return Δ_2d(n, m)
end

function Δ_1d_schurcompl(n)
    A = Δ_1d(2 * n)
    S = A[1:n, 1:n] - A[1:n, n+1:end] * inv(Matrix(A[n+1:end, n+1:end])) * A[n+1:end, 1:n]
    return S
end

function semisepbandedmatrix(
    n;
    r_lower = 1,
    r_upper = 1,
    band_lower = 1,
    band_upper = 1,
    seed = 1234,
)
    # set seed
    Random.seed!(seed)
    # intialize
    A = zeros(n, n)
    # semi separable part
    A =
        A +
        Array(LowerTriangular(rand(n, r_lower) * rand(r_lower, n))) +
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
    θ = collect(range(2 * pi / n, length = n, stop = 2 * pi))
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
