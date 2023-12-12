using Test
using Random
using LinearAlgebra


using SSSmatrices


include("casestudies.jl")

@testset "SSS constructor" begin
    #SSS
    no_blocks = 5
    n = [5, 5, 4, 4, 4]
    ranks_l = [2, 3, 2, 1]
    ranks_u = [3, 1, 2, 2]
    N = sum(n)
    D = [rand(n[i], n[i]) for i = 1:no_blocks]

    Q = [
        [rand(n[i], ranks_l[i]) for i = 1:no_blocks-1]
        [rand(n[no_blocks], 0)]
    ]
    R = [
        [rand(ranks_l[1], 0)]
        [rand(ranks_l[i+1], ranks_l[i]) for i = 1:no_blocks-2]
        [rand(0, ranks_l[no_blocks-1])]
    ]
    P = [
        [rand(n[1], 0)]
        [rand(n[i], ranks_l[i-1]) for i = 2:no_blocks]
    ]

    U = [
        [rand(n[i], ranks_u[i]) for i = 1:no_blocks-1]
        [rand(n[no_blocks], 0)]
    ]
    W = [
        [rand(ranks_u[1], 0)]
        [rand(ranks_u[i+1], ranks_u[i]) for i = 1:no_blocks-2]
        [rand(0, ranks_u[no_blocks-1])]
    ]
    V = [
        [rand(n[1], 0)]
        [rand(n[i], ranks_u[i-1]) for i = 2:no_blocks]
    ]

    A_SSS = SSS{Float64}(D, Q, R, P, U, W, V)
end

@testset "SSS matrix addition" begin
    A_SSS = random_SSS([5, 5, 4, 4, 4], [2, 3, 2, 1], [3, 1, 2, 2])
    B_SSS = random_SSS([5, 5, 4, 4, 4], [2, 3, 2, 1], [3, 1, 2, 2])
    @test Matrix(B_SSS + A_SSS) ≈ Matrix(B_SSS) + Matrix(A_SSS)
end

@testset "SSS matrix vector multiply" begin
    A_SSS = random_SSS([5, 5, 4, 4, 4], [2, 3, 2, 1], [3, 1, 2, 2])
    x = rand(size(A_SSS, 1))
    X = rand(size(A_SSS, 2), size(A_SSS, 2))
    Y = rand(15, size(A_SSS, 1))
    @test A_SSS * x ≈ Matrix(A_SSS) * x
    @test A_SSS * X ≈ Matrix(A_SSS) * X
end

@testset "Dense to SSS matrix construction" begin

    N = 200
    n = [10, 20, 40, 10, 60, 40, 20]
    tol = 1E-14

    # test on randomly generated matrix
    A = rand(N, N)
    A_SSS = SSS{Float64}(A, n, threshold=1E-13)
    rel_err = opnorm(Matrix(A_SSS) - A) / opnorm(A)
    @test rel_err < tol

    # test on 1D Laplacian
    A = Δ_1d(N)
    A_SSS = SSS{Float64}(A, n, threshold=1E-6)
    rel_err = opnorm(Matrix(A_SSS) - A) / opnorm(A)
    @test rel_err < tol
    #@test all(x -> x == 1, A_SSS.lower_hankel_ranks) & all(x -> x == 1, A_SSS.upper_hankel_ranks)

    # test on schur complement of 1D Laplacian 
    A = Δ_1d_schurcompl(N)
    A_SSS = SSS{Float64}(A, n, threshold=1E-8)
    rel_err = opnorm(Matrix(A_SSS) - A) / opnorm(A)
    @test rel_err < tol
    #@test all(x -> x == 1, A_SSS.lower_hankel_ranks) & all(x -> x == 1, A_SSS.upper_hankel_ranks)

    # test on semi separable banded matrix
    A = semisepbandedmatrix(N)
    A_SSS = SSS{Float64}(A, n, threshold=1E-8)
    rel_err = opnorm(Matrix(A_SSS) - A) / opnorm(A)
    @test rel_err < tol
    #@test all(x -> x == 1, A_SSS.lower_hankel_ranks) & all(x -> x == 2, A_SSS.upper_hankel_ranks)

    # # test construction Fourier Cauchy matrix
    # A = FourierCauchy{ComplexF64}(N)
    # A_SSS = SSS{ComplexF64}(A, determine_blocksizes_cauchy(N; K=3.0); threshold=1E-12)
    # rel_err = opnorm(Matrix(A_SSS) - A) / opnorm(A)
    # @test rel_err < tol

end

@testset "SSS solver" begin
    # 1 by 1 block SSS
    N = 30
    A = rand(N, N)
    x = rand(N)
    A_SSS = SSS{Float64}(A, [N], threshold = 1E-13)
    @test A \ x ≈ A_SSS \ x

    # 2 by 2 block SSS
    no_blocks = 2
    n = [5, 4]
    ranks_l = [2]
    ranks_u = [1]
    A_SSS = random_SSS(n, ranks_l, ranks_u)
    A = Matrix(A_SSS)
    c = rand(A_SSS.N)
    @test A \ c ≈ A_SSS \ c

    # 3 by 3 block SSS
    no_blocks = 3
    n = [5, 5, 4]
    ranks_l = [2, 2]
    ranks_u = [1, 2]
    A_SSS = random_SSS(n, ranks_l, ranks_u)
    A = Matrix(A_SSS)
    c = rand(A_SSS.N)
    @test A \ c ≈ A_SSS \ c

    # 4 by 4 block SSS
    no_blocks = 4
    n = [5, 3, 4, 5]
    ranks_l = [2, 1, 3]
    ranks_u = [1, 3, 2]
    A_SSS = random_SSS(n, ranks_l, ranks_u)
    A = Matrix(A_SSS)
    c = rand(A_SSS.N)
    @test A \ c ≈ A_SSS \ c

    # 8 by 8 block SSS
    no_blocks = 8
    n = [5, 3, 4, 5, 5, 5, 5, 5]
    ranks_l = [2, 1, 3, 2, 2, 1, 2]
    ranks_u = [1, 1, 2, 3, 2, 1, 2]
    A_SSS = random_SSS(n, ranks_l, ranks_u)
    A = Matrix(A_SSS)
    c = rand(A_SSS.N)
    @test A \ c ≈ A_SSS \ c

    # 10 by 10 block SSS
    no_blocks = 10
    n = [5, 3, 4, 5, 5, 5, 5, 5, 7, 6]
    ranks_l = [2, 1, 3, 2, 2, 1, 2, 3, 4]
    ranks_u = [1, 1, 2, 3, 2, 1, 2, 3, 4]
    A_SSS = random_SSS(n, ranks_l, ranks_u)
    A = Matrix(A_SSS)
    c = rand(A_SSS.N)
    @test A \ c ≈ A_SSS \ c

    N = 200
    n = [10, 20, 40, 10, 60, 40, 20]
    A = Δ_1d_schurcompl(N)
    A_SSS = SSS{Float64}(A, n)
    c = rand(N)

    @test A \ c ≈ A_SSS \ c

    # # FourierCauchy
    # N = 200
    # A = FourierCauchy{ComplexF64}(N)
    # A_SSS = SSS{ComplexF64}(A, determine_blocksizes_cauchy(N; K = 3.0); threshold = 1E-12)
    # c = rand(N)

    # @test Matrix(A_SSS) ≈ A
    # @test A \ c ≈ A_SSS \ c
end
