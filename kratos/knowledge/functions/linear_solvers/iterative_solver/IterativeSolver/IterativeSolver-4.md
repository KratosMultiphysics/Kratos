---
type: C++ Method
title: IterativeSolver
resource: linear_solvers/iterative_solver.h#L122-L126
generated:
  by: okf-rs/0.7.0
---

# Signature

`IterativeSolver(double NewTolerance, unsigned int NewMaxIterationsNumber, typename TPreconditionerType::Pointer pNewPreconditioner) : mResidualNorm(0), mIterationsNumber(0), mBNorm(0), mpPreconditioner(pNewPreconditioner), mTolerance(NewTolerance), mMaxIterationsNumber(NewMaxIterationsNumber)`