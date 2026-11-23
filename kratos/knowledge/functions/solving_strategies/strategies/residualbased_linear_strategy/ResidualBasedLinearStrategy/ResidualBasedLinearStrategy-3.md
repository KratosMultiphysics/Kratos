---
type: C++ Method
title: ResidualBasedLinearStrategy
resource: solving_strategies/strategies/residualbased_linear_strategy.h#L145-L185
generated:
  by: okf-rs/0.7.0
---

# Signature

`explicit ResidualBasedLinearStrategy( ModelPart& rModelPart, typename TSchemeType::Pointer pScheme, typename TLinearSolver::Pointer pNewLinearSolver, bool CalculateReactionFlag = false, bool ReformDofSetAtEachStep = false, bool CalculateNormDxFlag = false, bool MoveMeshFlag = false ) : BaseType(rModelPart, MoveMeshFlag), mpScheme(pScheme), mReformDofSetAtEachStep(ReformDofSetAtEachStep), mCalculateNormDxFlag(CalculateNormDxFlag), mCalculateReactionsFlag(CalculateReactionFlag)`