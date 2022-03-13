function tridiagonalexample(n)

    # x'' + 2x = f on [0,1], x(0) = x(1), x'(0) = x'(1)

    dx = 1.0 / n
    maindiag = (-2 / dx^2 + 2) * ones(n)
    subdiag = 1 / dx^2 * ones(n - 1)

    A = Array(Tridiagonal(subdiag, maindiag, subdiag))

    return A

end

function tridiagionalinverse(n)

    # x'' + 2x = f on [0,1], x(0) = x(1), x'(0) = x'(1)

    A = tridiagonalexample(n)
    A = inv(A)

    return A

end


function sanity_check1(n)

    # x'' + 2x = f on [0,1], x(0) = x(1), x'(0) = x'(1)

    dx = 1.0 / n
    maindiag = (-2 / dx^2 + 2) * ones(n)
    subdiag = 1 / dx^2 * ones(n - 1)

    A = Array(Tridiagonal(subdiag, maindiag, subdiag))
    A[1, n] = 1 / dx^2

    A[n, n] = -1 / (2 * dx)
    A[n, 2] = 1 / (2 * dx)
    A[n, n] = 0

    return A

end


function sanity_check2(n)

    sqrtn = Integer(ceil(sqrt(n)))

    twos = ones(n)
    minusones = -ones(n - 1)

    A = Array(Tridiagonal(minusones, twos, minusones))

    A[n+1-sqrtn:n, 1:sqrtn] += Array(I, sqrtn, sqrtn)
    A[1:sqrtn, n+1-sqrtn:n] += Array(I, sqrtn, sqrtn)

    return A

end

function sanity_check3(n)

    sqrtn = Integer(ceil(sqrt(n)))

    twos = ones(n)
    minusones = -ones(n - 1)

    A = Array(Tridiagonal(minusones, twos, minusones))
    A += 0.05 * rand(n, 1) * rand(1, n)

    A[n+1-sqrtn:n, 1:sqrtn] += Array(I, sqrtn, sqrtn)
    A[1:sqrtn, n+1-sqrtn:n] += Array(I, sqrtn, sqrtn)

    return A

end

function schur_comp_example(n)

    # x'' + 2x = f on [0,1], x(0) = x(1), x'(0) = x'(1)
    m = 2 * n

    dx = 1.0 / m
    maindiag = (-2 / dx^2 + 2) * ones(m)
    subdiag = 1 / dx^2 * ones(m - 1)

    A = Array(Tridiagonal(subdiag, maindiag, subdiag))
    A[1, m] = 1 / dx^2

    A[m, m] = -1 / (2 * dx)
    A[m, 2] = 1 / (2 * dx)
    A[m, m] = 0

    S = A[1:n, 1:n] - A[1:n, n+1:end] * inv(A[n+1:end, n+1:end]) * A[n+1:end, 1:n]

    return S

end

function schur_comp_example(n)

    # x'' + 2x = f on [0,1], x(0) = x(1), x'(0) = x'(1)
    m = 2 * n

    dx = 1.0 / m
    maindiag = (-2 / dx^2 + 2) * ones(m)
    subdiag = 1 / dx^2 * ones(m - 1)

    A = Array(Tridiagonal(subdiag, maindiag, subdiag))
    A[1, m] = 1 / dx^2

    A[m, m] = -1 / (2 * dx)
    A[m, 2] = 1 / (2 * dx)
    A[m, m] = 0

    println(inv(A[n+1:end, n+1:end]))

    S = A[1:n, 1:n] - A[1:n, n+1:end] * inv(A[n+1:end, n+1:end]) * A[n+1:end, 1:n]

    return S

end


function perturbedsemisepbandedmatrix(n; r_lower = 1, r_upper = 1, band_lower = 1,
    band_upper = 1, delta = 1, seed = 1234, corner_size = nothing)
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

    # Corner perturbations
    if corner_size == nothing
        corner_size = Integer(ceil(sqrt(n)))
    end

    A[n+1-corner_size:n, 1:corner_size] = A[n+1-corner_size:n, 1:corner_size] +
                                          delta * rand(corner_size, corner_size)
    A[1:corner_size, n+1-corner_size:n] = A[1:corner_size, n+1-corner_size:n] +
                                          delta * rand(corner_size, corner_size)

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


function circularlogkernel(n)

    # intialize
    A = zeros(n, n)

    for l = 0:n-1

        for k = 0:n-1

            if min(abs(k - l), abs(n - k + l)) == 1
                A[l+1, k+1] = (1 / n) * (1.5 + 1 / log(n)) * log((1 / n) * min(abs(k - l), abs(n - k + l)))
            elseif min(abs(k - l), abs(n - k + l)) > 1
                A[l+1, k+1] = (1 / n) * log((1 / n) * min(abs(k - l), abs(n - k + l)))
            end
        end

    end

    return A

end
