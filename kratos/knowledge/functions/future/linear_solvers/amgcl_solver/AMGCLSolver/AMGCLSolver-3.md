---
type: C++ Method
title: AMGCLSolver
resource: future/linear_solvers/amgcl_solver.h#L283-L310
generated:
  by: okf-rs/0.7.0
---

# Signature

`AMGCLSolver( AMGCLSmoother Smoother, AMGCLIterativeSolverType Solver, AMGCLCoarseningType Coarsening, double Tolerance, int MaxIterationsNumber, int Verbosity, int GMRESSize = 50, bool ProvideCoordinates = false ) : mTolerance(Tolerance), mMaxIterationsNumber(MaxIterationsNumber), mVerbosity(Verbosity), mBlockSize(1), mGMRESSize(GMRESSize), mCoarseEnough(1000), mFallbackToGMRES(false), mProvideCoordinates(ProvideCoordinates), mUseBlockMatricesIfPossible(false)`