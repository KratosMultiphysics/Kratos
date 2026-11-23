---
type: C++ Function
title: HasDofFor
resource: includes/node.h#L1367-L1375
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python_interface/kratos_globals/KratosGlobalsImpl/GetVariable
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/processes/apply_constant_vectorvalue_process/ApplyConstantVectorValueProcess/InternalApplyValue
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/flux_corrected_transport_convection_process/FluxCorrectedTransportConvectionProcess/ExecuteInitialize
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python/add_node_to_python/NodeHasDofFor
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/schemes/residual_based_bdf_custom_scheme/ResidualBasedBDFCustomScheme/InitializeSolutionStep
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/schemes/residual_based_bdf_custom_scheme/ComputePredictComponent
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/schemes/residual_based_bdf_displacement_scheme/ResidualBasedBDFDisplacementScheme/InitializeSolutionStep
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/schemes/residual_based_bdf_displacement_scheme/ResidualBasedBDFDisplacementScheme/Predict
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/schemes/residual_based_bossak_displacement_scheme/ResidualBasedBossakDisplacementScheme/Predict
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/schemes/residual_based_bossak_displacement_scheme/ResidualBasedBossakDisplacementScheme/InitializeSolutionStep
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/schemes/residual_based_pseudo_static_displacement_scheme/ResidualBasedPseudoStaticDisplacementScheme/Predict
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part/ModelPart/CreateNewMasterSlaveConstraint
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadConstraintsBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tensor_adaptors/fixity_tensor_adaptor/FixityTensorAdaptor/CheckContainer
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_specifications_utilities/TestSpecificationsUtilities/test_specifications_utilities_elements
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_specifications_utilities/TestSpecificationsUtilities/test_specifications_utilities_elements_list
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_specifications_utilities/TestSpecificationsUtilities/test_specifications_utilities_conditions
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_specifications_utilities/TestSpecificationsUtilities/test_specifications_utilities_conditions_list
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_specifications_utilities/TestSpecificationsUtilities/test_specifications_utilities_elements_dependencies
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_specifications_utilities/TestSpecificationsUtilities/test_specifications_utilities_elements_core_dependencies_list
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_specifications_utilities/TestSpecificationsUtilities/test_specifications_utilities_elements_dependencies_list
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_specifications_utilities/TestSpecificationsUtilities/test_specifications_utilities_conditions_dependencies
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_specifications_utilities/TestSpecificationsUtilities/test_specifications_utilities_conditions_dependencies_list
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_variable_utils/TestVariableUtils/testAddDofsList
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_variable_utils/TestVariableUtils/testAddDofsWithReactionList
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/assign_master_slave_constraints_to_neighbours_utility/AssignMasterSlaveConstraintsToNeighboursUtility/GetDofsAndCoordinatesForSlaveNode
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/assign_master_slave_constraints_to_neighbours_utility/AssignMasterSlaveConstraintsToNeighboursUtility/GetDofsAndCoordinatesForCloudOfNodes
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/specifications_utilities/SpecificationsUtilities/AddMissingDofsFromSpecifications
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/variable_utils/ApplyFixity
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`inline bool HasDofFor(const VariableData& rDofVariable) const`

# Calls

- [GetVariable](../../../functions/python_interface/kratos_globals/KratosGlobalsImpl/GetVariable.md)

# Called by

- [InternalApplyValue](../../../functions/processes/apply_constant_vectorvalue_process/ApplyConstantVectorValueProcess/InternalApplyValue.md)
- [ExecuteInitialize](../../../functions/processes/flux_corrected_transport_convection_process/FluxCorrectedTransportConvectionProcess/ExecuteInitialize.md)
- [NodeHasDofFor](../../../functions/python/add_node_to_python/NodeHasDofFor.md)
- [InitializeSolutionStep](../../../functions/solving_strategies/schemes/residual_based_bdf_custom_scheme/ResidualBasedBDFCustomScheme/InitializeSolutionStep.md)
- [ComputePredictComponent](../../../functions/solving_strategies/schemes/residual_based_bdf_custom_scheme/ComputePredictComponent.md)
- [InitializeSolutionStep](../../../functions/solving_strategies/schemes/residual_based_bdf_displacement_scheme/ResidualBasedBDFDisplacementScheme/InitializeSolutionStep.md)
- [Predict](../../../functions/solving_strategies/schemes/residual_based_bdf_displacement_scheme/ResidualBasedBDFDisplacementScheme/Predict.md)
- [Predict](../../../functions/solving_strategies/schemes/residual_based_bossak_displacement_scheme/ResidualBasedBossakDisplacementScheme/Predict.md)
- [InitializeSolutionStep](../../../functions/solving_strategies/schemes/residual_based_bossak_displacement_scheme/ResidualBasedBossakDisplacementScheme/InitializeSolutionStep.md)
- [Predict](../../../functions/solving_strategies/schemes/residual_based_pseudo_static_displacement_scheme/ResidualBasedPseudoStaticDisplacementScheme/Predict.md)
- [CreateNewMasterSlaveConstraint](../../../functions/sources/model_part/ModelPart/CreateNewMasterSlaveConstraint.md)
- [ReadConstraintsBlock](../../../functions/sources/model_part_io/ModelPartIO/ReadConstraintsBlock.md)
- [CheckContainer](../../../functions/tensor_adaptors/fixity_tensor_adaptor/FixityTensorAdaptor/CheckContainer.md)
- [test_specifications_utilities_elements](../../../functions/tests/test_specifications_utilities/TestSpecificationsUtilities/test_specifications_utilities_elements.md)
- [test_specifications_utilities_elements_list](../../../functions/tests/test_specifications_utilities/TestSpecificationsUtilities/test_specifications_utilities_elements_list.md)
- [test_specifications_utilities_conditions](../../../functions/tests/test_specifications_utilities/TestSpecificationsUtilities/test_specifications_utilities_conditions.md)
- [test_specifications_utilities_conditions_list](../../../functions/tests/test_specifications_utilities/TestSpecificationsUtilities/test_specifications_utilities_conditions_list.md)
- [test_specifications_utilities_elements_dependencies](../../../functions/tests/test_specifications_utilities/TestSpecificationsUtilities/test_specifications_utilities_elements_dependencies.md)
- [test_specifications_utilities_elements_core_dependencies_list](../../../functions/tests/test_specifications_utilities/TestSpecificationsUtilities/test_specifications_utilities_elements_core_dependencies_list.md)
- [test_specifications_utilities_elements_dependencies_list](../../../functions/tests/test_specifications_utilities/TestSpecificationsUtilities/test_specifications_utilities_elements_dependencies_list.md)
- [test_specifications_utilities_conditions_dependencies](../../../functions/tests/test_specifications_utilities/TestSpecificationsUtilities/test_specifications_utilities_conditions_dependencies.md)
- [test_specifications_utilities_conditions_dependencies_list](../../../functions/tests/test_specifications_utilities/TestSpecificationsUtilities/test_specifications_utilities_conditions_dependencies_list.md)
- [testAddDofsList](../../../functions/tests/test_variable_utils/TestVariableUtils/testAddDofsList.md)
- [testAddDofsWithReactionList](../../../functions/tests/test_variable_utils/TestVariableUtils/testAddDofsWithReactionList.md)
- [GetDofsAndCoordinatesForSlaveNode](../../../functions/utilities/assign_master_slave_constraints_to_neighbours_utility/AssignMasterSlaveConstraintsToNeighboursUtility/GetDofsAndCoordinatesForSlaveNode.md)
- [GetDofsAndCoordinatesForCloudOfNodes](../../../functions/utilities/assign_master_slave_constraints_to_neighbours_utility/AssignMasterSlaveConstraintsToNeighboursUtility/GetDofsAndCoordinatesForCloudOfNodes.md)
- [AddMissingDofsFromSpecifications](../../../functions/utilities/specifications_utilities/SpecificationsUtilities/AddMissingDofsFromSpecifications.md)
- [ApplyFixity](../../../functions/utilities/variable_utils/ApplyFixity.md)