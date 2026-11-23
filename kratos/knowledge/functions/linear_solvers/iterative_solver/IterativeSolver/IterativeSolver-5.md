---
type: C++ Method
title: IterativeSolver
resource: linear_solvers/iterative_solver.h#L128-L155
generated:
  by: okf-rs/0.7.0
---

# Signature

`IterativeSolver(Parameters settings, typename TPreconditionerType::Pointer pNewPreconditioner = Kratos::make_shared<TPreconditionerType>() ): mResidualNorm(0), mIterationsNumber(0), mBNorm(0), mpPreconditioner(pNewPreconditioner)`