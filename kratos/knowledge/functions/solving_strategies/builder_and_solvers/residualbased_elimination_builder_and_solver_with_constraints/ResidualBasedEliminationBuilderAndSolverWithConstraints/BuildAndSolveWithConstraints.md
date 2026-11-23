---
type: C++ Method
title: BuildAndSolveWithConstraints
resource: solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_with_constraints.h#L458-L516
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_with_constraints/ResidualBasedEliminationBuilderAndSolverWithConstraints/ApplyMasterSlaveRelation
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_with_constraints/ResidualBasedEliminationBuilderAndSolverWithConstraints/ComputeEffectiveConstant
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_with_constraints/ResidualBasedEliminationBuilderAndSolverWithConstraints/BuildWithConstraints
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/timer/Timer/Stop
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/builtin_timer/BuiltinTimer/BuiltinTimer
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/timer/Timer/Start
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_with_constraints/ResidualBasedEliminationBuilderAndSolverWithConstraints/ReconstructSlaveSolutionAfterSolve
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/communicator/Communicator/MyPID
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_with_constraints/ResidualBasedEliminationBuilderAndSolverWithConstraints/BuildAndSolve
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void BuildAndSolveWithConstraints( typename TSchemeType::Pointer pScheme, ModelPart& rModelPart, TSystemMatrixType& rA, TSystemVectorType& rDx, TSystemVectorType& rb )`

# Calls

- [ApplyMasterSlaveRelation](../../../../../functions/solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_with_constraints/ResidualBasedEliminationBuilderAndSolverWithConstraints/ApplyMasterSlaveRelation.md)
- [ComputeEffectiveConstant](../../../../../functions/solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_with_constraints/ResidualBasedEliminationBuilderAndSolverWithConstraints/ComputeEffectiveConstant.md)
- [BuildWithConstraints](../../../../../functions/solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_with_constraints/ResidualBasedEliminationBuilderAndSolverWithConstraints/BuildWithConstraints.md)
- [Stop](../../../../../functions/utilities/timer/Timer/Stop.md)
- [BuiltinTimer](../../../../../functions/utilities/builtin_timer/BuiltinTimer/BuiltinTimer.md)
- [Start](../../../../../functions/utilities/timer/Timer/Start.md)
- [ReconstructSlaveSolutionAfterSolve](../../../../../functions/solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_with_constraints/ResidualBasedEliminationBuilderAndSolverWithConstraints/ReconstructSlaveSolutionAfterSolve.md)
- [MyPID](../../../../../functions/sources/communicator/Communicator/MyPID.md)

# Called by

- [BuildAndSolve](../../../../../functions/solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_with_constraints/ResidualBasedEliminationBuilderAndSolverWithConstraints/BuildAndSolve.md)