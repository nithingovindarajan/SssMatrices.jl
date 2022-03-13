using Test
using Random
using LinearAlgebra


using SSSmatrices


greet()


@testset "SSS matrices sanity check" begin


    # Construct SSS matrix directly from generators
    N = 20
    n = [5, 5, 4, 4, 4, 4, 6, 7, 8, 6, 5, 5, 4, 4, 4, 4, 6, 7, 8, 6]
    Gpi = [2, 3, 2, 2, 2, 2, 2, 2, 2, 2, 3, 2, 2, 2, 2, 2, 2, 2, 3]
    Hpi = [4, 2, 2, 4, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2, 2, 2, 2, 2, 4]


    Di = [rand(n[i], n[i]) for i = 1:N]
    Ui = [rand(n[i], Hpi[i]) for i = 1:N-1]
    Wi = [rand(Hpi[i], Hpi[i+1]) for i = 1:N-2]
    for i = 1:N-2
        Wi[i] = Wi[i] / (2 * norm(Wi[i]))
    end
    Vi = [rand(n[i+1], Hpi[i]) for i = 1:N-1]
    Pi = [rand(n[i+1], Gpi[i]) for i = 1:N-1]
    Ri = [rand(Gpi[i+1], Gpi[i]) for i = 1:N-2]
    for i = 1:N-2
        Ri[i] = Ri[i] / (2 * norm(Ri[i]))
    end
    Qi = [rand(n[i], Gpi[i]) for i = 1:N-1]

    A = SSS{Float64}(N, n, Gpi, Hpi, Di, Ui, Wi, Vi, Pi, Ri, Qi)

    Adiag = DiagonalSSS{Float64}(N, n, Gpi, Hpi, Di, Ui, Wi, Vi, Pi, Ri, Qi)

    # SSS matrix construction from dense matrix
    n = [5, 5, 5, 5, 5, 5, 5, 5]
    N = sum(n)
    A = rand(N, N)
    A_SSS = SSS(A, n)
    Adense = Matrix(A_SSS)

    @test Adense ≈ A


    # SSS multiply
    b = rand(sum(n))
    c_SSS = A_SSS * b
    c = A * b

    @test c_SSS ≈ c              # Error  original vs SSS


    # SSS solve
    x_SSS = A_SSS \ b
    x = A \ b

    @test x_SSS ≈ x              # Error SSS solver vs Dense solver

end;








