# SssMatrices.jl

`SssMatrices` is Julia package for Sequentially Semi-Separable (SSS) matrices. SSS matrices were originally introduced in [1] within the context systems theory. This package implements routines to convert a dense matrix into (minimal) SSS representations and to solve linear systems Ax = b in SSS form. 

## Getting started

Directly construct SSS matrices from its generators

```Julia

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

# construct SSS matrix
A_SSS = SSS{Float64}(D, Q, R, P, U, W, V)
```

Performing algebra: adding two SSS matrices yields another SSS matrix, right multiplication of an SSS matrix
with a vector calls up specialized routines.

```Julia

# randomly generate an SSS matrix of a specific rank profile.
A_SSS = random_SSS([5, 5, 4, 4, 4], [2, 3, 2, 1], [3, 1, 2, 2])
B_SSS = random_SSS([5, 5, 4, 4, 4], [2, 3, 2, 1], [3, 1, 2, 2])

# adding two SSS matrix
C_SSS = B_SSS + A_SSS
Matrix(C_SSS) ≈ Matrix(B_SSS) + Matrix(A_SSS) # evaluates to true

# Matrix vector multiplication
x = rand(size(A_SSS, 1))
y = A_SSS * x
y ≈ Matrix(A_SSS) * x # evaluates to true
``````

We may converting a dense matrix into SSS form and solve it with specialized solvers.

```Julia
using LinearAlgebra
using BlockDiagonals
using SssMatrices

# dimensions
n = [20, 20, 30, 10, 20]
N = sum(n)

# dense matrix
D = BlockDiagonal([rand(k, k) for k in n])
L = LowerTriangular(rand(N,1) * rand(1,N))
U = UpperTriangular(rand(N,1) * rand(1,N))
A = D + L + U

# convert to SSS form
A_SSS = SSS{Float64}(A, n, threshold=1E-13)

# solve
b = rand(N)
x_SSS = A_SSS \ b
x = A \ b
x ≈ x_SSS # evaluates to true
```

## Relevant literature

[1] Dewilde, Patrick, and Alle-Jan Van der Veen. Time-varying systems and computations. Springer Science & Business Media, 1998.

[2] Chandrasekaran, S., Dewilde, P., Gu, M., Pals, T., Sun, X., van der Veen, A. J., & White, D. (2005). Some fast algorithms for sequentially semiseparable representations. SIAM Journal on Matrix Analysis and Applications, 27(2), 341-364.



<!-- 
GOALS:
- QR-based solver
- add ULV
- Transposing SSS matrices
- LU decomposition
- QR decomposition
- SSS matrix-matrix multiply
- SSS inverse
- recompress operation
- extnd solvers to minimum norm and overdetermined systems
- Toeplitz solver: Chandrasekaran, S., Gu, M., Sun, X., Xia, J., & Zhu, J. (2008). A superfast algorithm for Toeplitz systems of linear equations. SIAM Journal on Matrix Analysis and Applications, 29(4), 1247-1266.
 -->





