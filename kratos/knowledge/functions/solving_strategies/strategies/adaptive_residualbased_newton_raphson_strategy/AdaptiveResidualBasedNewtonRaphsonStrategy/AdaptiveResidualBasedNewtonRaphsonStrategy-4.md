---
type: C++ Method
title: AdaptiveResidualBasedNewtonRaphsonStrategy
resource: solving_strategies/strategies/adaptive_residualbased_newton_raphson_strategy.h#L257-L319
generated:
  by: okf-rs/0.7.0
---

# Signature

`explicit AdaptiveResidualBasedNewtonRaphsonStrategy( ModelPart& rModelPart, typename TSchemeType::Pointer pScheme, typename TLinearSolver::Pointer pNewLinearSolver, typename TConvergenceCriteriaType::Pointer pNewConvergenceCriteria, typename TBuilderAndSolverType::Pointer pNewBuilderAndSolver, int MaxIterations = 30, int MinIterations = 4, bool CalculateReactions = false, bool ReformDofSetAtEachStep = false, bool MoveMeshFlag = false, double ReductionFactor = 0.5, double IncreaseFactor = 1.3, int NumberOfCycles = 5 ) : BaseType(rModelPart, MoveMeshFlag)`