---
type: C++ Method
title: ImplicitStrategy
resource: future/solving_strategies/strategies/implicit_strategy.h#L132-L152
generated:
  by: okf-rs/0.7.0
---

# Signature

`explicit ImplicitStrategy( ModelPart& rModelPart, SchemePointerType pScheme, LinearSolverPointerType pLinearSolver, bool ComputeReactions = false, bool ReformDofSetAtEachStep = false, bool CalculateNormDxFlag = false, bool MoveMeshFlag = false) : Strategy(rModelPart, Parameters()) , mpScheme(pScheme) , mpLinearSolver(pLinearSolver) , mReformDofsAtEachStep(ReformDofSetAtEachStep) , mComputeReactions(ComputeReactions)`