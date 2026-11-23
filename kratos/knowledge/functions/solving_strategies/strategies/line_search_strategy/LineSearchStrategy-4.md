---
type: C++ Function
title: LineSearchStrategy
resource: solving_strategies/strategies/line_search_strategy.h#L167-L181
generated:
  by: okf-rs/0.7.0
---

# Signature

`explicit LineSearchStrategy( ModelPart& rModelPart, typename TSchemeType::Pointer pScheme, typename TConvergenceCriteriaType::Pointer pNewConvergenceCriteria, typename TBuilderAndSolverType::Pointer pNewBuilderAndSolver, int MaxIterations = 30, bool CalculateReactions = false, bool ReformDofSetAtEachStep = false, bool MoveMeshFlag = false ) : BaseType(rModelPart, pScheme, pNewConvergenceCriteria, pNewBuilderAndSolver, MaxIterations, CalculateReactions, ReformDofSetAtEachStep, MoveMeshFlag)`