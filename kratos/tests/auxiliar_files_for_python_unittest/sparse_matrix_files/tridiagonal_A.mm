%%MatrixMarket matrix coordinate real general
%=================================================================================
%
% A small tridiagonal system, of the shape a 1D transient groundwater flow problem
% on three two-node line elements produces: a permeability (1D Laplacian) operator
% scaled by k/mu = 9.084e-06 / 1.0e-03, with the first degree of freedom prescribed
% (its row zeroed off the diagonal and its column zeroed in the other rows, as
% ApplyDirichletConditions leaves it - the zeroed entries are kept explicitly so the
% sparsity pattern matches the assembled one).
%
% The point of this matrix is that it is tridiagonal, so ILU0 drops no fill-in and
% is an *exact* factorization of it. A Krylov method preconditioned with it reaches
% the solution within the first iteration, which is the degenerate case ("happy"
% breakdown) that iterative solvers have to commit rather than discard.
%
%=================================================================================
4 4 10
1 1 9.084e-03
1 2 0.0
2 1 0.0
2 2 1.8168e-02
2 3 -9.084e-03
3 2 -9.084e-03
3 3 1.8168e-02
3 4 -9.084e-03
4 3 -9.084e-03
4 4 9.084e-03
