---
type: C++ Method
title: AdjointStrategy
resource: tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme.cpp#L104-L119
generated:
  by: okf-rs/0.7.0
---

# Signature

`AdjointStrategy(ModelPart& rModelPart, Kratos::shared_ptr<PrimalResults> pPrimalResults, AdjointResponseFunction::Pointer pResponseFunction) : ImplicitSolvingStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType>(rModelPart), mpPrimalResults(pPrimalResults)`