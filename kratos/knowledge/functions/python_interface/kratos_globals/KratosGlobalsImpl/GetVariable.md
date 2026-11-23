---
type: Python Method
title: GetVariable
resource: python_interface/kratos_globals.py#L35-L71
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/constraints/slip_constraint/SlipConstraint/PrintInfo
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/dof/Dof/SetNodalData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/dof/Dof/Info
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/dof/Dof/PrintData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/node/Fix
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/node/Free
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/node/GetDofPosition
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/node/pGetDof
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/node/pAddDof
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/node/HasDofFor
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/node/IsFixed
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/node/SortDofs
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/linear_solvers/amgcl_ns_solver/AMGCL_NS_Solver/ProvideAdditionalData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/structured_mesh_refinement_modeler/GenerateNode
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/assign_average_master_slave_constraints_process/AssignAverageMasterSlaveConstraintsProcess/ExecuteInitialize
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/assign_master_slave_constraints_to_neighbours_process/AssignMasterSlaveConstraintsToNeighboursProcess/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/assign_scalar_variable_process/AssignScalarVariableProcess/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/assign_time_derivative_process/AssignTimeDerivativeProcess/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/assign_vector_variable_process/AssignVectorVariableProcess/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/assign_vector_variable_to_entities_process/AssignVectorVariableToEntitiesProcess/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/auxiliary_solver_utilities/AddVariables
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/auxiliary_solver_utilities/AddDofs
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/check_scalar_base_process/CheckScalarBaseProcess/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/fix_scalar_variable_process/FixScalarVariableProcess/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/fix_vector_variable_process/FixVectorVariableProcess/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/gid_hdf5_postprocess/GiDHDF5OutputProcess/_GetListOfVariables
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/gid_output/GiDOutput/write_results
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/json_output_process/JsonOutputProcess/__generate_variable_list_from_input
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/kratos_utilities/GenerateVariableListFromInput
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/point_output_process/PointOutputProcess/ExecuteInitialize
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/print_info_in_file_process/PrintInfoInFileProcess/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/vtu_output_process/VtuOutputProcess/ExecuteInitialize
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/p_multigrid/p_multigrid_builder_and_solver/PMultigridBuilderAndSolver<TSparse,TDense>/SetUpDofSet
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/p_multigrid/p_multigrid_utilities/MakePRestrictionOperator
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/p_multigrid/sparse_utilities/ApplyDirichletConditions
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/p_multigrid/status_stream/MakeVtuOutput
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/convergencecriterias/mixed_generic_criteria/MixedGenericCriteria/FindVarLocalKey
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/strategies/residualbased_newton_raphson_strategy/WriteDofInfo
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/WriteConstraints
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/sources/test_dof/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/sources/test_node/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/strategies/builder_and_solvers/test_explicit_builder_and_solver/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/utilities/test_dof_array_utilities/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_constraint_restart/TestConstraintRestart/test_constraint_serialization
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_dofs/TestDofs/testDofExport
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_kratos_globals/TestKratosGlobals/test_get_variable
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/constraint_utilities/PreComputeExplicitConstraintConstribution
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/constraint_utilities/PreComputeExplicitConstraintMassAndInertia
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def GetVariable(self, VarName):`

# Called by

- [PrintInfo](../../../../functions/constraints/slip_constraint/SlipConstraint/PrintInfo.md)
- [SetNodalData](../../../../functions/includes/dof/Dof/SetNodalData.md)
- [Info](../../../../functions/includes/dof/Dof/Info.md)
- [PrintData](../../../../functions/includes/dof/Dof/PrintData.md)
- [Fix](../../../../functions/includes/node/Fix.md)
- [Free](../../../../functions/includes/node/Free.md)
- [GetDofPosition](../../../../functions/includes/node/GetDofPosition.md)
- [pGetDof](../../../../functions/includes/node/pGetDof.md)
- [pAddDof](../../../../functions/includes/node/pAddDof.md)
- [HasDofFor](../../../../functions/includes/node/HasDofFor.md)
- [IsFixed](../../../../functions/includes/node/IsFixed.md)
- [SortDofs](../../../../functions/includes/node/SortDofs.md)
- [ProvideAdditionalData](../../../../functions/linear_solvers/amgcl_ns_solver/AMGCL_NS_Solver/ProvideAdditionalData.md)
- [GenerateNode](../../../../functions/modeler/structured_mesh_refinement_modeler/GenerateNode.md)
- [ExecuteInitialize](../../../../functions/python_scripts/assign_average_master_slave_constraints_process/AssignAverageMasterSlaveConstraintsProcess/ExecuteInitialize.md)
- [__init__](../../../../functions/python_scripts/assign_master_slave_constraints_to_neighbours_process/AssignMasterSlaveConstraintsToNeighboursProcess/__init__.md)
- [__init__](../../../../functions/python_scripts/assign_scalar_variable_process/AssignScalarVariableProcess/__init__.md)
- [__init__](../../../../functions/python_scripts/assign_time_derivative_process/AssignTimeDerivativeProcess/__init__.md)
- [__init__](../../../../functions/python_scripts/assign_vector_variable_process/AssignVectorVariableProcess/__init__.md)
- [__init__](../../../../functions/python_scripts/assign_vector_variable_to_entities_process/AssignVectorVariableToEntitiesProcess/__init__.md)
- [AddVariables](../../../../functions/python_scripts/auxiliary_solver_utilities/AddVariables.md)
- [AddDofs](../../../../functions/python_scripts/auxiliary_solver_utilities/AddDofs.md)
- [__init__](../../../../functions/python_scripts/check_scalar_base_process/CheckScalarBaseProcess/__init__.md)
- [__init__](../../../../functions/python_scripts/fix_scalar_variable_process/FixScalarVariableProcess/__init__.md)
- [__init__](../../../../functions/python_scripts/fix_vector_variable_process/FixVectorVariableProcess/__init__.md)
- [_GetListOfVariables](../../../../functions/python_scripts/gid_hdf5_postprocess/GiDHDF5OutputProcess/_GetListOfVariables.md)
- [write_results](../../../../functions/python_scripts/gid_output/GiDOutput/write_results.md)
- [__generate_variable_list_from_input](../../../../functions/python_scripts/json_output_process/JsonOutputProcess/__generate_variable_list_from_input.md)
- [GenerateVariableListFromInput](../../../../functions/python_scripts/kratos_utilities/GenerateVariableListFromInput.md)
- [ExecuteInitialize](../../../../functions/python_scripts/point_output_process/PointOutputProcess/ExecuteInitialize.md)
- [__init__](../../../../functions/python_scripts/print_info_in_file_process/PrintInfoInFileProcess/__init__.md)
- [ExecuteInitialize](../../../../functions/python_scripts/vtu_output_process/VtuOutputProcess/ExecuteInitialize.md)
- [SetUpDofSet](../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/p_multigrid_builder_and_solver/PMultigridBuilderAndSolver<TSparse,TDense>/SetUpDofSet.md)
- [MakePRestrictionOperator](../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/p_multigrid_utilities/MakePRestrictionOperator.md)
- [ApplyDirichletConditions](../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/sparse_utilities/ApplyDirichletConditions.md)
- [MakeVtuOutput](../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/status_stream/MakeVtuOutput.md)
- [FindVarLocalKey](../../../../functions/solving_strategies/convergencecriterias/mixed_generic_criteria/MixedGenericCriteria/FindVarLocalKey.md)
- [WriteDofInfo](../../../../functions/solving_strategies/strategies/residualbased_newton_raphson_strategy/WriteDofInfo.md)
- [WriteConstraints](../../../../functions/sources/model_part_io/ModelPartIO/WriteConstraints.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/sources/test_dof/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/sources/test_node/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/strategies/builder_and_solvers/test_explicit_builder_and_solver/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/utilities/test_dof_array_utilities/KRATOS_TEST_CASE_IN_SUITE.md)
- [test_constraint_serialization](../../../../functions/tests/test_constraint_restart/TestConstraintRestart/test_constraint_serialization.md)
- [testDofExport](../../../../functions/tests/test_dofs/TestDofs/testDofExport.md)
- [test_get_variable](../../../../functions/tests/test_kratos_globals/TestKratosGlobals/test_get_variable.md)
- [PreComputeExplicitConstraintConstribution](../../../../functions/utilities/constraint_utilities/PreComputeExplicitConstraintConstribution.md)
- [PreComputeExplicitConstraintMassAndInertia](../../../../functions/utilities/constraint_utilities/PreComputeExplicitConstraintMassAndInertia.md)