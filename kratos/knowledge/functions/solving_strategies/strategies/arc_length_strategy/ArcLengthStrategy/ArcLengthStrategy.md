---
type: C++ Method
title: ArcLengthStrategy
resource: solving_strategies/strategies/arc_length_strategy.h#L89-L102
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/GetInt
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetBool
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`ArcLengthStrategy( ModelPart& model_part, typename TSchemeType::Pointer pScheme, typename TConvergenceCriteriaType::Pointer pNewConvergenceCriteria, typename TBuilderAndSolverType::Pointer pNewBuilderAndSolver, Parameters ThisParameters ) : ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(model_part, pScheme, pNewConvergenceCriteria, pNewBuilderAndSolver, ThisParameters["max_iteration"].GetInt(), ThisParameters["compute_reactions"].GetBool(), ThisParameters["reform_dofs_at_each_step"].GetBool(), ThisParameters["move_mesh_flag"].GetBool())`

# Calls

- [GetInt](../../../../../functions/sources/kratos_parameters/Parameters/GetInt.md)
- [GetBool](../../../../../functions/sources/kratos_parameters/Parameters/GetBool.md)