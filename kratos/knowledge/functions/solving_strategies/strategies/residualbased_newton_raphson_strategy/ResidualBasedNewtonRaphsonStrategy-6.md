---
type: C++ Function
title: ResidualBasedNewtonRaphsonStrategy
resource: solving_strategies/strategies/residualbased_newton_raphson_strategy.h#L278-L302
generated:
  by: okf-rs/0.7.0
---

# Signature

`explicit ResidualBasedNewtonRaphsonStrategy( ModelPart& rModelPart, typename TSchemeType::Pointer pScheme, typename TLinearSolver::Pointer pNewLinearSolver, typename TConvergenceCriteriaType::Pointer pNewConvergenceCriteria, typename TBuilderAndSolverType::Pointer pNewBuilderAndSolver, int MaxIterations = 30, bool CalculateReactions = false, bool ReformDofSetAtEachStep = false, bool MoveMeshFlag = false) : ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(rModelPart, pScheme, pNewConvergenceCriteria, pNewBuilderAndSolver, MaxIterations, CalculateReactions, ReformDofSetAtEachStep, MoveMeshFlag)`