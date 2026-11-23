---
type: C++ Function
title: ResidualBasedNewtonRaphsonStrategy
resource: solving_strategies/strategies/residualbased_newton_raphson_strategy.h#L312-L348
generated:
  by: okf-rs/0.7.0
---

# Signature

`ResidualBasedNewtonRaphsonStrategy( ModelPart& rModelPart, typename TSchemeType::Pointer pScheme, typename TLinearSolver::Pointer pNewLinearSolver, typename TConvergenceCriteriaType::Pointer pNewConvergenceCriteria, Parameters Settings) : BaseType(rModelPart), mpScheme(pScheme), mpConvergenceCriteria(pNewConvergenceCriteria), mInitializeWasPerformed(false), mKeepSystemConstantDuringIterations(false)`