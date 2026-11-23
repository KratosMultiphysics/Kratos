---
type: C++ Function
title: SolutionStepsDataHas
resource: includes/node.h#L631-L634
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/elements/distance_calculation_element_simplex/DistanceCalculationElementSimplex/Check
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/elements/edge_based_gradient_recovery_element/EdgeBasedGradientRecoveryElement<TDim>/Check
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/elements/embedded_nodal_variable_calculation_element_simplex/EmbeddedNodalVariableCalculationElementSimplex<double>/Check
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/elements/embedded_nodal_variable_calculation_element_simplex/EmbeddedNodalVariableCalculationElementSimplex<array_1d<double,3>>/Check
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/utilities/parallel_fill_communicator/ParallelFillCommunicator/ComputeCommunicationPlan
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_discontinuous_distance_to_skin_process/CalculateEmbeddedVariableFromSkinSpecialization
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_nodal_area_process/CalculateNodalAreaProcess<THistorical>/Execute
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/find_nodal_h_process/if
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/parallel_distance_calculation_process/ParallelDistanceCalculationProcess<TDim>/Check
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/tetrahedral_mesh_orientation_check/TetrahedralMeshOrientationCheck/Execute
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python/add_node_to_python/NodeSolutionStepsDataHas
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python/add_node_to_python/IndexingUtility
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/table_accessor/TableAccessor/GetValueFromTable
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/processes/test_calculate_distance_to_skin_process/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_model_part/TestModelPart/test_modelpart_variables_list
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_time_averaging/TimeAveragingProcessTests/__CreateModel
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/coordinate_transformation_utilities/CalculateRotationOperatorPureShapeSensitivities
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/embedded_skin_utility/InterpolateMeshVariableToSkinSpecialization
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/search_utilities/PreparePointsSearch
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/variable_utils/AddDof
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/variable_utils/AddDofWithReaction
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/variable_utils/HasValue
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`bool SolutionStepsDataHas(const VariableData& rThisVariable) const`

# Called by

- [Check](../../../functions/elements/distance_calculation_element_simplex/DistanceCalculationElementSimplex/Check.md)
- [Check](../../../functions/elements/edge_based_gradient_recovery_element/EdgeBasedGradientRecoveryElement<TDim>/Check.md)
- [Check](../../../functions/elements/embedded_nodal_variable_calculation_element_simplex/EmbeddedNodalVariableCalculationElementSimplex<double>/Check.md)
- [Check](../../../functions/elements/embedded_nodal_variable_calculation_element_simplex/EmbeddedNodalVariableCalculationElementSimplex<array_1d<double,3>>/Check.md)
- [ComputeCommunicationPlan](../../../functions/mpi/utilities/parallel_fill_communicator/ParallelFillCommunicator/ComputeCommunicationPlan.md)
- [CalculateEmbeddedVariableFromSkinSpecialization](../../../functions/processes/calculate_discontinuous_distance_to_skin_process/CalculateEmbeddedVariableFromSkinSpecialization.md)
- [Execute](../../../functions/processes/calculate_nodal_area_process/CalculateNodalAreaProcess<THistorical>/Execute.md)
- [if](../../../functions/processes/find_nodal_h_process/if.md)
- [Check](../../../functions/processes/parallel_distance_calculation_process/ParallelDistanceCalculationProcess<TDim>/Check.md)
- [Execute](../../../functions/processes/tetrahedral_mesh_orientation_check/TetrahedralMeshOrientationCheck/Execute.md)
- [NodeSolutionStepsDataHas](../../../functions/python/add_node_to_python/NodeSolutionStepsDataHas.md)
- [IndexingUtility](../../../functions/python/add_node_to_python/IndexingUtility.md)
- [GetValueFromTable](../../../functions/sources/table_accessor/TableAccessor/GetValueFromTable.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../functions/tests/cpp_tests/processes/test_calculate_distance_to_skin_process/KRATOS_TEST_CASE_IN_SUITE.md)
- [test_modelpart_variables_list](../../../functions/tests/test_model_part/TestModelPart/test_modelpart_variables_list.md)
- [__CreateModel](../../../functions/tests/test_time_averaging/TimeAveragingProcessTests/__CreateModel.md)
- [CalculateRotationOperatorPureShapeSensitivities](../../../functions/utilities/coordinate_transformation_utilities/CalculateRotationOperatorPureShapeSensitivities.md)
- [InterpolateMeshVariableToSkinSpecialization](../../../functions/utilities/embedded_skin_utility/InterpolateMeshVariableToSkinSpecialization.md)
- [PreparePointsSearch](../../../functions/utilities/search_utilities/PreparePointsSearch.md)
- [AddDof](../../../functions/utilities/variable_utils/AddDof.md)
- [AddDofWithReaction](../../../functions/utilities/variable_utils/AddDofWithReaction.md)
- [HasValue](../../../functions/utilities/variable_utils/HasValue.md)