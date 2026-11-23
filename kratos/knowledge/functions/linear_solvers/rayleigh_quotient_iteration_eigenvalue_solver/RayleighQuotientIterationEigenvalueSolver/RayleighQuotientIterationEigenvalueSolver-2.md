---
type: C++ Method
title: RayleighQuotientIterationEigenvalueSolver
resource: linear_solvers/rayleigh_quotient_iteration_eigenvalue_solver.h#L111-L123
generated:
  by: okf-rs/0.7.0
---

# Signature

`RayleighQuotientIterationEigenvalueSolver( double NewMaxTolerance, unsigned int NewMaxIterationsNumber, unsigned int NewRequiredEigenvalueNumber, typename TLinearSolverType::Pointer pLinearSolver, double ShiftingConvergence = 0.25 ): BaseType(NewMaxTolerance, NewMaxIterationsNumber), mRequiredEigenvalueNumber(NewRequiredEigenvalueNumber), mpLinearSolver(pLinearSolver), mShiftingConvergence(ShiftingConvergence)`