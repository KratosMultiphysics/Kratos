---
type: C++ Method
title: PMGStatusStream
resource: solving_strategies/builder_and_solvers/p_multigrid/status_stream.cpp#L286-L302
visibility: private
generated:
  by: okf-rs/0.7.0
---

# Signature

`PMGStatusStream::PMGStatusStream(std::ostream& rStream, PMultigridBuilderAndSolver<TSparse,TDense>& rBuilderAndSolver, ModelPart& rModelPart, const TMatrix& rRootLhs, const TVector& rRootSolution, const TVector& rRootRhs, bool UseAnsiColors) : mpImpl(new Impl {/*mpStream=*/ &rStream, /*mpModelPart=*/ &rModelPart, /*mSystem=*/ std::make_tuple(&rRootLhs, &rRootSolution, &rRootRhs), /*mpBuilderAndSolver=*/ &rBuilderAndSolver, /*mpMaybeVtuOutput=*/ {}, /*mVtuPaths=*/ {}, /*mMaybeLastReport=*/ {}, /*mUseAnsiColors=*/ UseAnsiColors})`