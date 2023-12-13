# SssMatrices.jl

`SssMatrices` is Julia package for Sequentially Semi-Separable (SSS) matrices. SSS matrices were originally introduced in [1] within the context systems theory. This package implements routines to convert a dense matrix into (minimal) SSS representations and to solve linear systems Ax = b in SSS form. 

## Getting started


```Julia
using LinearAlgebra
using SssMatrices

K(x,y) = (x-y) != 0 ? 1/(x-y) : 1.
A = [ K(x,y) for x=-1:0.001:1, y=-1:0.001:1]
A_SSS = SSS(A)



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
- extnd solvers to minimum norm and overdetermined systems
- Toeplitz solver: Chandrasekaran, S., Gu, M., Sun, X., Xia, J., & Zhu, J. (2008). A superfast algorithm for Toeplitz systems of linear equations. SIAM Journal on Matrix Analysis and Applications, 29(4), 1247-1266.
 -->





