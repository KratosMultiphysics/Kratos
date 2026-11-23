---
type: C++ Method
title: ResidualBasedLinearStrategy
resource: solving_strategies/strategies/residualbased_linear_strategy.h#L198-L233
generated:
  by: okf-rs/0.7.0
---

# Signature

`explicit ResidualBasedLinearStrategy( ModelPart& rModelPart, typename TSchemeType::Pointer pScheme, typename TBuilderAndSolverType::Pointer pNewBuilderAndSolver, bool CalculateReactionFlag = false, bool ReformDofSetAtEachStep = false, bool CalculateNormDxFlag = false, bool MoveMeshFlag = false ) : BaseType(rModelPart, MoveMeshFlag), mpScheme(pScheme), mpBuilderAndSolver(pNewBuilderAndSolver), mReformDofSetAtEachStep(ReformDofSetAtEachStep), mCalculateNormDxFlag(CalculateNormDxFlag), mCalculateReactionsFlag(CalculateReactionFlag)`