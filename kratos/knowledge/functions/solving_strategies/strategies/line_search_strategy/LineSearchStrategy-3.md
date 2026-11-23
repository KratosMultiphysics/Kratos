---
type: C++ Function
title: LineSearchStrategy
resource: solving_strategies/strategies/line_search_strategy.h#L140-L154
generated:
  by: okf-rs/0.7.0
---

# Signature

`explicit LineSearchStrategy( ModelPart& rModelPart, typename TSchemeType::Pointer pScheme, typename TLinearSolver::Pointer pNewLinearSolver, typename TConvergenceCriteriaType::Pointer pNewConvergenceCriteria, int MaxIterations = 30, bool CalculateReactions = false, bool ReformDofSetAtEachStep = false, bool MoveMeshFlag = false ) : BaseType(rModelPart, pScheme, pNewLinearSolver,pNewConvergenceCriteria,MaxIterations,CalculateReactions,ReformDofSetAtEachStep, MoveMeshFlag)`