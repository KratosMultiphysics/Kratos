---
type: C++ Method
title: AddValue
resource: sources/kratos_parameters.cpp#L455-L466
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/input_output/cad_json_output/CadJsonOutput/GetParameters
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/cad_json_output/CadJsonOutput/GetBrepSurfaceParameters
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/cad_json_output/CadJsonOutput/GetBoundaryLoopParameters
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/serial_model_part_combinator_modeler/SerialModelPartCombinatorModeler/SetupModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/sub_model_part_skin_detection_process/SubModelPartSkinDetectionProcess<TDim>/GetDefaultParameters
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/apply_periodic_condition_process/ApplyPeriodicConditionProcessPythonWrapper/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/assign_mpcs_to_neighbours_process/AssignMPCsToNeighboursProcess/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/assign_scalar_input_to_entities_process/AssignScalarInputToEntitiesProcess/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/assign_scalar_variable_to_entities_process/AssignScalarVariableToEntitiesProcess/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/assign_vector_by_direction_process/AssignVectorByDirectionProcess/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/assign_vector_by_direction_to_entity_process/AssignVectorByDirectionToEntityProcess/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/assign_vector_variable_process/AssignVectorVariableProcess/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/assign_vector_variable_to_entities_process/AssignVectorVariableToEntitiesProcess/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/check_scalar_to_nodes_process/CheckScalarToNodesProcess/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/check_vector_components_to_nodes_process/CheckVectorComponentsToNodesProcess/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/fix_vector_variable_process/FixVectorVariableProcess/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/gid_hdf5_postprocess/GiDHDF5OutputProcess/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/gid_output_process/GiDOutputProcess/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/integration_values_extrapolation_to_nodes_process/IntegrationValuesExtrapolationToNodesProcess/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/print_info_in_file_process/PrintInfoInFileProcess/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/python_solver/PythonSolver/_GetRestartSettings
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/restart_utility/RestartUtility/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/save_restart_process/SaveRestartProcess/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/skin_detection_process/SkinDetectionProcess/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/time_based_ascii_file_writer_utility/TimeBasedAsciiFileWriterUtility/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/timer_process/Factory
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/p_multigrid/augmented_lagrange_constraint_assembler/AugmentedLagrangeConstraintAssembler<TSparse,TDense>/AugmentedLagrangeConstraintAssembler
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/p_multigrid/master_slave_constraint_assembler/MasterSlaveConstraintAssembler<TSparse,TDense>/MasterSlaveConstraintAssembler
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/p_multigrid/p_multigrid_builder_and_solver/PMultigridBuilderAndSolver<TSparse,TDense>/AssignSettings
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/AddDouble
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/AddInt
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/AddBool
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/AddString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/AddStringArray
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/AddVector
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/AddMatrix
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/CopyValuesFromExistingParameters
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/assign_unique_model_part_collection_tag_utility/AssignUniqueModelPartCollectionTagUtility/WriteTagsToJson
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/interval_utility/IntervalUtility/IntervalUtility
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/profiler/Profiler/Write
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/sensitivity_builder/ReplaceDeprecatedNameIfExists
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/specifications_utilities/SpecificationsUtilities/GetDocumention
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void Parameters::AddValue( const std::string& rEntry, const Parameters& rOtherValue )`

# Called by

- [GetParameters](../../../../functions/input_output/cad_json_output/CadJsonOutput/GetParameters.md)
- [GetBrepSurfaceParameters](../../../../functions/input_output/cad_json_output/CadJsonOutput/GetBrepSurfaceParameters.md)
- [GetBoundaryLoopParameters](../../../../functions/input_output/cad_json_output/CadJsonOutput/GetBoundaryLoopParameters.md)
- [SetupModelPart](../../../../functions/modeler/serial_model_part_combinator_modeler/SerialModelPartCombinatorModeler/SetupModelPart.md)
- [GetDefaultParameters](../../../../functions/processes/sub_model_part_skin_detection_process/SubModelPartSkinDetectionProcess<TDim>/GetDefaultParameters.md)
- [__init__](../../../../functions/python_scripts/apply_periodic_condition_process/ApplyPeriodicConditionProcessPythonWrapper/__init__.md)
- [__init__](../../../../functions/python_scripts/assign_mpcs_to_neighbours_process/AssignMPCsToNeighboursProcess/__init__.md)
- [__init__](../../../../functions/python_scripts/assign_scalar_input_to_entities_process/AssignScalarInputToEntitiesProcess/__init__.md)
- [__init__](../../../../functions/python_scripts/assign_scalar_variable_to_entities_process/AssignScalarVariableToEntitiesProcess/__init__.md)
- [__init__](../../../../functions/python_scripts/assign_vector_by_direction_process/AssignVectorByDirectionProcess/__init__.md)
- [__init__](../../../../functions/python_scripts/assign_vector_by_direction_to_entity_process/AssignVectorByDirectionToEntityProcess/__init__.md)
- [__init__](../../../../functions/python_scripts/assign_vector_variable_process/AssignVectorVariableProcess/__init__.md)
- [__init__](../../../../functions/python_scripts/assign_vector_variable_to_entities_process/AssignVectorVariableToEntitiesProcess/__init__.md)
- [__init__](../../../../functions/python_scripts/check_scalar_to_nodes_process/CheckScalarToNodesProcess/__init__.md)
- [__init__](../../../../functions/python_scripts/check_vector_components_to_nodes_process/CheckVectorComponentsToNodesProcess/__init__.md)
- [__init__](../../../../functions/python_scripts/fix_vector_variable_process/FixVectorVariableProcess/__init__.md)
- [__init__](../../../../functions/python_scripts/gid_hdf5_postprocess/GiDHDF5OutputProcess/__init__.md)
- [__init__](../../../../functions/python_scripts/gid_output_process/GiDOutputProcess/__init__.md)
- [__init__](../../../../functions/python_scripts/integration_values_extrapolation_to_nodes_process/IntegrationValuesExtrapolationToNodesProcess/__init__.md)
- [__init__](../../../../functions/python_scripts/print_info_in_file_process/PrintInfoInFileProcess/__init__.md)
- [_GetRestartSettings](../../../../functions/python_scripts/python_solver/PythonSolver/_GetRestartSettings.md)
- [__init__](../../../../functions/python_scripts/restart_utility/RestartUtility/__init__.md)
- [__init__](../../../../functions/python_scripts/save_restart_process/SaveRestartProcess/__init__.md)
- [__init__](../../../../functions/python_scripts/skin_detection_process/SkinDetectionProcess/__init__.md)
- [__init__](../../../../functions/python_scripts/time_based_ascii_file_writer_utility/TimeBasedAsciiFileWriterUtility/__init__.md)
- [Factory](../../../../functions/python_scripts/timer_process/Factory.md)
- [AugmentedLagrangeConstraintAssembler](../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/augmented_lagrange_constraint_assembler/AugmentedLagrangeConstraintAssembler<TSparse,TDense>/AugmentedLagrangeConstraintAssembler.md)
- [MasterSlaveConstraintAssembler](../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/master_slave_constraint_assembler/MasterSlaveConstraintAssembler<TSparse,TDense>/MasterSlaveConstraintAssembler.md)
- [AssignSettings](../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/p_multigrid_builder_and_solver/PMultigridBuilderAndSolver<TSparse,TDense>/AssignSettings.md)
- [AddDouble](../../../../functions/sources/kratos_parameters/Parameters/AddDouble.md)
- [AddInt](../../../../functions/sources/kratos_parameters/Parameters/AddInt.md)
- [AddBool](../../../../functions/sources/kratos_parameters/Parameters/AddBool.md)
- [AddString](../../../../functions/sources/kratos_parameters/Parameters/AddString.md)
- [AddStringArray](../../../../functions/sources/kratos_parameters/Parameters/AddStringArray.md)
- [AddVector](../../../../functions/sources/kratos_parameters/Parameters/AddVector.md)
- [AddMatrix](../../../../functions/sources/kratos_parameters/Parameters/AddMatrix.md)
- [CopyValuesFromExistingParameters](../../../../functions/sources/kratos_parameters/Parameters/CopyValuesFromExistingParameters.md)
- [WriteTagsToJson](../../../../functions/utilities/assign_unique_model_part_collection_tag_utility/AssignUniqueModelPartCollectionTagUtility/WriteTagsToJson.md)
- [IntervalUtility](../../../../functions/utilities/interval_utility/IntervalUtility/IntervalUtility.md)
- [Write](../../../../functions/utilities/profiler/Profiler/Write.md)
- [ReplaceDeprecatedNameIfExists](../../../../functions/utilities/sensitivity_builder/ReplaceDeprecatedNameIfExists.md)
- [GetDocumention](../../../../functions/utilities/specifications_utilities/SpecificationsUtilities/GetDocumention.md)