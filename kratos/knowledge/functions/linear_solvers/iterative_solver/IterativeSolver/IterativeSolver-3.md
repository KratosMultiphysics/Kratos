---
type: C++ Method
title: IterativeSolver
resource: linear_solvers/iterative_solver.h#L114-L120
generated:
  by: okf-rs/0.7.0
---

# Signature

`IterativeSolver(double NewTolerance, unsigned int NewMaxIterationsNumber) : mResidualNorm(0) , mIterationsNumber(0) , mBNorm(0) , mpPreconditioner(new TPreconditionerType()) , mTolerance(NewTolerance) , mMaxIterationsNumber(NewMaxIterationsNumber)`