---
type: C++ Function
title: ResidualBasedNewtonRaphsonStrategy
resource: solving_strategies/strategies/residualbased_newton_raphson_strategy.h#L359-L397
generated:
  by: okf-rs/0.7.0
---

# Signature

`ResidualBasedNewtonRaphsonStrategy( ModelPart& rModelPart, typename TSchemeType::Pointer pScheme, typename TConvergenceCriteriaType::Pointer pNewConvergenceCriteria, typename TBuilderAndSolverType::Pointer pNewBuilderAndSolver, Parameters Settings) : BaseType(rModelPart), mpScheme(pScheme), mpBuilderAndSolver(pNewBuilderAndSolver), mpConvergenceCriteria(pNewConvergenceCriteria), mInitializeWasPerformed(false), mKeepSystemConstantDuringIterations(false)`