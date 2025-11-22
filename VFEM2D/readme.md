About the matlab libarary Verified Finite Element Method (VFEM2D) 

VFEM2D aims to provide rigorous matrices assembling in the FEM computation.
Currently, this libarary provide the matrices computation for the following FEMs.
- Langrange FEM of arbitrary order
- Raviart-Thomas FEM space of arbitrary order.

Also, it provides the algorithm to obtain rigorous eigenvalue bound computation for the Laplace operator and the Steklov operator. Liu's projection-based eigenvalue estimation algorithm [3,4] along with linear conforming FEM is utilized to obtain gurantteed lower eigenvalue bound. Also, the Lehmann-Goerisch method (Chap. 5 of [1]) is used to obtain high-precision eigenvaue bounds along with the Lagrange FEM of higher degrees.

Revision history:
- This mathlab library was originally developed by Chun'guang You and Xuefeng Liu in 2016 for the rigorous eigenvalue estimation for Stekolv eigenvalue bounds [2].
- 2024: Xuefeng LIU, revision for Dirichlet eigenvalue computation.
- 2025/11/21: Xuefeng LIu, sorted version for release.

Reference
1. Xuefeng LIU, Guaranteed Computational Methods for Self-Adjoint Differential Eigenvalue Problems, 2024, SpringerBriefs in Mathematics, Springer Singapore.
2. Chun'guang You, Hehu Xie and Xuefeng Liu, Guaranteed Eigenvalue Bounds for the Steklov Eigenvalue Problem, SIAM Journal on Numerical Analysis, 57(3), 1395-1410, 2019. DOI:10.1137/18M1189592.
3. Xuefeng Liu, A framework of verified eigenvalue bounds for self-adjoint differential operators, Applied Mathematics and Computation, 267, pp.341-355, 2015. DOI:10.1016/j.amc.2015.03.048.
4. Xuefeng Liu, Shin'ichi Oishi, Verified eigenvalue evaluation for the Laplacian over polygonal domains of arbitrary shape, SIAM Journal on Numerical Analysis, 51(3), 1634-1654, 2013. DOI:10.1137/120878446. 