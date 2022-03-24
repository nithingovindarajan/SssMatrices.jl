using Test
using Random
using LinearAlgebra


using SSSmatrices


greet()


@testset "SSStypes constructors" begin


    # Diagonal SSS
    no_blocks = 5
    n = [5, 5, 4, 4, 4]
    N = sum(n)
    D = [rand(n[i], n[i]) for i = 1:no_blocks]

    D_SSS = DiagonalSSS(D)


    @test D_SSS.no_blocks == no_blocks
    @test D_SSS.n == n
    @test D_SSS.N == N
    @test D_SSS.diagonal.D == D


    #SSS
    no_blocks = 5
    n = [5, 5, 4, 4, 4]
    ranks_l = [2, 3, 2, 1]
    ranks_u = [3, 1, 2, 2]
    N = sum(n)
    D = [rand(n[i], n[i]) for i = 1:no_blocks]

    Q = [[rand(n[i], ranks_l[i]) for i = 1:no_blocks-1]
        [rand(n[no_blocks], 0)]]
    R = [[rand(ranks_l[1], 0)]
        [rand(ranks_l[i+1], ranks_l[i]) for i = 1:no_blocks-2]
        [rand(0, ranks_l[no_blocks-1])]]
    P = [[rand(n[1], 0)]
        [rand(n[i], ranks_l[i-1]) for i = 2:no_blocks]]

    U = [[rand(n[i], ranks_u[i]) for i = 1:no_blocks-1]
        [rand(n[no_blocks], 0)]]
    W = [[rand(ranks_u[1], 0)]
        [rand(ranks_u[i+1], ranks_u[i]) for i = 1:no_blocks-2]
        [rand(0, ranks_u[no_blocks-1])]]
    V = [[rand(n[1], 0)]
        [rand(n[i], ranks_u[i-1]) for i = 2:no_blocks]]

    A_SSS = SSS(D, Q, R, P, U, W, V)

    x = rand(size(A_SSS, 1))
    X = rand(size(A_SSS, 2), 20)
    Y = rand(15, size(A_SSS, 1))

    @test A_SSS * x ≈ Matrix(A_SSS) * x
    @test A_SSS * X ≈ Matrix(A_SSS) * X

    # TODO: @test Y * A_SSS ≈ Y * Matrix(A_SSS)


end

@testset "SSS matrix vector multiply" begin


    A_SSS = random_SSS([5, 5, 4, 4, 4], [2, 3, 2, 1], [3, 1, 2, 2])

    x = rand(size(A_SSS, 1))
    X = rand(size(A_SSS, 2), size(A_SSS, 2))
    Y = rand(15, size(A_SSS, 1))

    @test A_SSS * x ≈ Matrix(A_SSS) * x
    @test A_SSS * X ≈ Matrix(A_SSS) * X
    @test A_SSS * X' ≈ Matrix(A_SSS) * X'
    @test A_SSS' * X ≈ Matrix(A_SSS)' * X
    @test Y * A_SSS ≈ Y * Matrix(A_SSS)
    @test Y * A_SSS' ≈ Y * Matrix(A_SSS)'


end



@testset "SSS matrix addition" begin

    A_SSS = random_SSS([5, 5, 4, 4, 4], [2, 3, 2, 1], [3, 1, 2, 2])
    B_SSS = random_SSS([5, 5, 4, 4, 4], [2, 3, 2, 1], [3, 1, 2, 2])

    @test Matrix(B_SSS + A_SSS) ≈ Matrix(B_SSS) + Matrix(A_SSS)
    @test Matrix(B_SSS' + A_SSS) ≈ Matrix(B_SSS)' + Matrix(A_SSS)
    @test Matrix(B_SSS + A_SSS') ≈ Matrix(B_SSS) + Matrix(A_SSS)'
    @test Matrix(B_SSS' + A_SSS') ≈ Matrix(B_SSS)' + Matrix(A_SSS)'
end



@testset "Fast LDR matrix solvers" begin

    n = 100
    b = rand(n)
    A = SquareToeplitz(rand(n), rand(n - 1))

    @test Matrix(A) \ b ≈ A \ b

end



@testset "Dense to SSS matrix construction" begin


    # test on randomly generated matrix
    A = rand(200, 200)
    n = [10, 20, 40, 10, 60, 40, 20]

    A_SSS = SSS(A, n)
    A_dense = Matrix(A)

    @test A_dense ≈ A



end


