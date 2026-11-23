---
type: C++ Function
title: LineSearchStrategy
resource: solving_strategies/strategies/line_search_strategy.h#L196-L211
generated:
  by: okf-rs/0.7.0
---

# Signature

`explicit LineSearchStrategy( ModelPart& rModelPart, typename TSchemeType::Pointer pScheme, typename TLinearSolver::Pointer pNewLinearSolver, typename TConvergenceCriteriaType::Pointer pNewConvergenceCriteria, typename TBuilderAndSolverType::Pointer pNewBuilderAndSolver, int MaxIterations = 30, bool CalculateReactions = false, bool ReformDofSetAtEachStep = false, bool MoveMeshFlag = false ) : BaseType(rModelPart, pScheme, pNewLinearSolver,pNewConvergenceCriteria,pNewBuilderAndSolver,MaxIterations,CalculateReactions,ReformDofSetAtEachStep, MoveMeshFlag)`