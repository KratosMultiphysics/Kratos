---
type: C++ Method
title: MyPID
resource: sources/communicator.cpp#L91-L94
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/input_output/ensight_output/EnSightOutput/GetOutputFileName
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/vtk_output/VtkOutput/GetOutputFileName
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/vtu_output/VtuOutput/PrintOutput
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/tests/cpp_tests/sources/test_mpi_communicator/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/tests/test_mpi_communicator/TestMPICommunicator/testCommunicatorRankSize
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/tests/test_mpi_model_part/TestMPIModelPart/test_remove_nodes_parallel_interfaces
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/ensight_output_process/EnSightOutputProcess/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/tikz_output_process/TikZOutputProcess/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/vtk_output_process/VtkOutputProcess/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/explicit_builder/ExplicitBuilder/SetUpDofSet
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver/ResidualBasedEliminationBuilderAndSolver/SystemSolveWithPhysics
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_with_constraints/ResidualBasedEliminationBuilderAndSolverWithConstraints/BuildAndSolveWithConstraints
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_with_constraints/ResidualBasedEliminationBuilderAndSolverWithConstraints/SetUpDofSetWithConstraints
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_with_constraints/ResidualBasedEliminationBuilderAndSolverWithConstraints/SystemSolveWithPhysics
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/strategies/residualbased_newton_raphson_strategy/EchoInfo
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/dof_utilities/block_build_dof_array_utility/BlockBuildDofArrayUtility/SetUpDofArray
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/dof_utilities/dof_array_utilities/DofArrayUtilities/SetUpDofArray
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/dof_utilities/dof_array_utilities/DofArrayUtilities/SetUpEffectiveDofArray
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/dof_utilities/elimination_build_dof_array_utility/EliminationBuildDofArrayUtility/SetUpDofArray
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`int Communicator::MyPID() const`

# Called by

- [GetOutputFileName](../../../../functions/input_output/ensight_output/EnSightOutput/GetOutputFileName.md)
- [GetOutputFileName](../../../../functions/input_output/vtk_output/VtkOutput/GetOutputFileName.md)
- [PrintOutput](../../../../functions/input_output/vtu_output/VtuOutput/PrintOutput.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/mpi/tests/cpp_tests/sources/test_mpi_communicator/KRATOS_TEST_CASE_IN_SUITE.md)
- [testCommunicatorRankSize](../../../../functions/mpi/tests/test_mpi_communicator/TestMPICommunicator/testCommunicatorRankSize.md)
- [test_remove_nodes_parallel_interfaces](../../../../functions/mpi/tests/test_mpi_model_part/TestMPIModelPart/test_remove_nodes_parallel_interfaces.md)
- [__init__](../../../../functions/python_scripts/ensight_output_process/EnSightOutputProcess/__init__.md)
- [__init__](../../../../functions/python_scripts/tikz_output_process/TikZOutputProcess/__init__.md)
- [__init__](../../../../functions/python_scripts/vtk_output_process/VtkOutputProcess/__init__.md)
- [SetUpDofSet](../../../../functions/solving_strategies/builder_and_solvers/explicit_builder/ExplicitBuilder/SetUpDofSet.md)
- [SystemSolveWithPhysics](../../../../functions/solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver/ResidualBasedEliminationBuilderAndSolver/SystemSolveWithPhysics.md)
- [BuildAndSolveWithConstraints](../../../../functions/solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_with_constraints/ResidualBasedEliminationBuilderAndSolverWithConstraints/BuildAndSolveWithConstraints.md)
- [SetUpDofSetWithConstraints](../../../../functions/solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_with_constraints/ResidualBasedEliminationBuilderAndSolverWithConstraints/SetUpDofSetWithConstraints.md)
- [SystemSolveWithPhysics](../../../../functions/solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_with_constraints/ResidualBasedEliminationBuilderAndSolverWithConstraints/SystemSolveWithPhysics.md)
- [EchoInfo](../../../../functions/solving_strategies/strategies/residualbased_newton_raphson_strategy/EchoInfo.md)
- [SetUpDofArray](../../../../functions/utilities/dof_utilities/block_build_dof_array_utility/BlockBuildDofArrayUtility/SetUpDofArray.md)
- [SetUpDofArray](../../../../functions/utilities/dof_utilities/dof_array_utilities/DofArrayUtilities/SetUpDofArray.md)
- [SetUpEffectiveDofArray](../../../../functions/utilities/dof_utilities/dof_array_utilities/DofArrayUtilities/SetUpEffectiveDofArray.md)
- [SetUpDofArray](../../../../functions/utilities/dof_utilities/elimination_build_dof_array_utility/EliminationBuildDofArrayUtility/SetUpDofArray.md)