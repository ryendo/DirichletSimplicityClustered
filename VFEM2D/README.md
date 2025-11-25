Verified Finite Element Method (2D)
-------------------------------

VFEM2D is a MATLAB library for rigorous finite element matrix assembly and verified eigenvalue computation on planar domains. It focuses on producing mathematically guaranteed results suitable for verified numerics and certified error bounds.

Key features
- High-order Lagrange finite elements (arbitrary polynomial degree)
- Raviart–Thomas finite elements (arbitrary order)
- Rigorous assembly routines designed for verified linear algebra
- Verified eigenvalue bounds for Laplace and Steklov operators

Algorithms and methods
- Lower eigenvalue bounds: projection-based verified estimators (Liu’s framework [1,3,4])
- High-precision bounds: Lehmann–Goerisch method combined with higher-order conforming FEM (Chap 5 of [1])
- Support for both Dirichlet eigenvalue problems and Steklov eigenvalue problems

Revision history
- Original MATLAB implementation (2016) — Chun'guang You and Xuefeng Liu (Steklov eigenvalue verification [2])
- 2024 — Xuefeng Liu: revisions for Dirichlet eigenvalue computations
- 2025-11-21 — Xuefeng Liu: sorted and packaged release

References
1. X. Liu, Guaranteed Computational Methods for Self-Adjoint Differential Eigenvalue Problems, SpringerBriefs in Mathematics, 2024. (https://doi.org/10.1007/978-981-97-3577-8) 
2. C. You, H. Xie, X. Liu, Guaranteed Eigenvalue Bounds for the Steklov Eigenvalue Problem, SIAM J. Numer. Anal., 2019. (DOI:10.1137/18M1189592) 
3. X. Liu, A framework of verified eigenvalue bounds for self-adjoint differential operators, Appl. Math. Comput., 2015.  (DOI:10.1016/j.amc.2015.03.048)
4. X. Liu, S. Oishi, Verified eigenvalue evaluation for the Laplacian over polygonal domains, SIAM J. Numer. Anal., 2013.　(DOI:10.1137/120878446)

Third-party software
- VEIGS (verified eigenvalue solver for sparse matrices, MATLAB): included submodule (version committed 2025-11-21, commit 99dd50a). Repository: https://github.com/yuuka-math/veigs