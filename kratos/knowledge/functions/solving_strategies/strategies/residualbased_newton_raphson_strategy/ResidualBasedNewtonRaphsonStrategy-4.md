---
type: C++ Function
title: ResidualBasedNewtonRaphsonStrategy
resource: solving_strategies/strategies/residualbased_newton_raphson_strategy.h#L166-L208
generated:
  by: okf-rs/0.7.0
---

# Signature

`explicit ResidualBasedNewtonRaphsonStrategy( ModelPart& rModelPart, typename TSchemeType::Pointer pScheme, typename TLinearSolver::Pointer pNewLinearSolver, typename TConvergenceCriteriaType::Pointer pNewConvergenceCriteria, int MaxIterations = 30, bool CalculateReactions = false, bool ReformDofSetAtEachStep = false, bool MoveMeshFlag = false) : BaseType(rModelPart, MoveMeshFlag), mpScheme(pScheme), mpConvergenceCriteria(pNewConvergenceCriteria), mReformDofSetAtEachStep(ReformDofSetAtEachStep), mCalculateReactionsFlag(CalculateReactions), mMaxIterationNumber(MaxIterations), mInitializeWasPerformed(false), mKeepSystemConstantDuringIterations(false)`