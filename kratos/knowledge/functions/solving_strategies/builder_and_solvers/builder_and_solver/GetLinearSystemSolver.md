---
type: C++ Function
title: GetLinearSystemSolver
resource: solving_strategies/builder_and_solvers/builder_and_solver.h#L247-L250
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/solving_strategies/builder_and_solvers/p_multigrid/p_multigrid_builder_and_solver/PMultigridBuilderAndSolver<TSparse,TDense>/Impl/ExecuteMultigridLoop
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/p_multigrid/p_multigrid_builder_and_solver/PMultigridBuilderAndSolver<TSparse,TDense>/Impl/ExecuteConstraintLoop
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/strategies/residualbased_newton_raphson_strategy/ResidualBasedNewtonRaphsonStrategy
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`typename TLinearSolver::Pointer GetLinearSystemSolver() const`

# Called by

- [ExecuteMultigridLoop](../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/p_multigrid_builder_and_solver/PMultigridBuilderAndSolver<TSparse,TDense>/Impl/ExecuteMultigridLoop.md)
- [ExecuteConstraintLoop](../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/p_multigrid_builder_and_solver/PMultigridBuilderAndSolver<TSparse,TDense>/Impl/ExecuteConstraintLoop.md)
- [ResidualBasedNewtonRaphsonStrategy](../../../../functions/solving_strategies/strategies/residualbased_newton_raphson_strategy/ResidualBasedNewtonRaphsonStrategy.md)