---
type: C++ Function
title: IsNot
resource: containers/flags.h#L335-L338
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/elements/distance_calculation_element_simplex/DistanceCalculationElementSimplex/CalculateLocalSystem
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/clean_up_problematic_triangles_modeler/CleanUpProblematicTrianglesModeler/CleanUpProblematicGeometries
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/clean_up_problematic_triangles_modeler/CleanUpProblematicTrianglesModeler/ComputeNullAreaTriangles
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_embedded_signed_distance_to_3d_skin_process/CalculateEmbeddedSignedDistanceTo3DSkinProcess/PeakValuesCorrection
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/MappingPressureToStructure
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/fast_transfer_between_model_parts_process/FastTransferBetweenModelPartsProcess/Execute
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/from_json_check_result_process/FromJSONCheckResultProcess/InitializeDatabases
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/from_json_check_result_process/FromJSONCheckResultProcess/SizeDatabase
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/from_json_check_result_process/CheckFlag
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/mesh_coarsening_process/MeshCoarseningProcess/SelectCoarseMeshNodes
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/mesh_local_smoothing_process/MeshLocalSmoothingProcess/PerformSmoothing
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/mesh_node_collapsing_process/MeshNodeCollapsingProcess/CollapseNode
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/parallel_distance_calculation_process/ParallelDistanceCalculationProcess<TDim>/SetDistancesOnDividedElements
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/parallel_distance_calculation_process/ParallelDistanceCalculationProcess<TDim>/ExtendDistancesByLayer
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/tetrahedral_mesh_orientation_check/TetrahedralMeshOrientationCheck/Execute
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver/ResidualBasedBlockBuilderAndSolver/InternalSystemSolveWithPhysics
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver/ResidualBasedBlockBuilderAndSolver/BuildMasterSlaveConstraints
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver_with_lagrange_multiplier/ResidualBasedBlockBuilderAndSolverWithLagrangeMultiplier/ApplyRHSConstraints
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part/RemoveEntities
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ModelPartIO
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/WriteModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadSubModelPartBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/processes/test_find_intersected_geometrical_objects_process/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_flags/TestFlags/testFlagSetReset
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_variable_utils/TestVariableUtils/test_set_flag
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/auxiliar_model_part_utilities/AuxiliarModelPartUtilities/RemoveElementAndBelongings
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/auxiliar_model_part_utilities/AuxiliarModelPartUtilities/RemoveElementsAndBelongings
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/auxiliar_model_part_utilities/AuxiliarModelPartUtilities/RemoveConditionAndBelongings
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/auxiliar_model_part_utilities/AuxiliarModelPartUtilities/RemoveConditionsAndBelongings
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/convect_particles_utilities/ParticleConvectUtily/EraseOuterElements
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/exact_mortar_segmentation_utility/ExactMortarIntegrationUtility<TDim, TNumNodes, TBelong, TNumNodesMaster>/TestIODebug
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/mortar_utilities/ComputeTangentNodeWithLMAndSlip
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/shifted_boundary_meshless_interface_utility/ShiftedBoundaryMeshlessInterfaceUtility/CalculateGradientBasedConformingExtensionBasis
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/shifted_boundary_meshless_interface_utility/ShiftedBoundaryMeshlessInterfaceUtility/CalculateMeshlessBasedConformingExtensionBasis
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`bool IsNot(Flags const& rOther) const`

# Called by

- [CalculateLocalSystem](../../../functions/elements/distance_calculation_element_simplex/DistanceCalculationElementSimplex/CalculateLocalSystem.md)
- [CleanUpProblematicGeometries](../../../functions/modeler/clean_up_problematic_triangles_modeler/CleanUpProblematicTrianglesModeler/CleanUpProblematicGeometries.md)
- [ComputeNullAreaTriangles](../../../functions/modeler/clean_up_problematic_triangles_modeler/CleanUpProblematicTrianglesModeler/ComputeNullAreaTriangles.md)
- [PeakValuesCorrection](../../../functions/processes/calculate_embedded_signed_distance_to_3d_skin_process/CalculateEmbeddedSignedDistanceTo3DSkinProcess/PeakValuesCorrection.md)
- [MappingPressureToStructure](../../../functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/MappingPressureToStructure.md)
- [Execute](../../../functions/processes/fast_transfer_between_model_parts_process/FastTransferBetweenModelPartsProcess/Execute.md)
- [InitializeDatabases](../../../functions/processes/from_json_check_result_process/FromJSONCheckResultProcess/InitializeDatabases.md)
- [SizeDatabase](../../../functions/processes/from_json_check_result_process/FromJSONCheckResultProcess/SizeDatabase.md)
- [CheckFlag](../../../functions/processes/from_json_check_result_process/CheckFlag.md)
- [SelectCoarseMeshNodes](../../../functions/processes/mesh_coarsening_process/MeshCoarseningProcess/SelectCoarseMeshNodes.md)
- [PerformSmoothing](../../../functions/processes/mesh_local_smoothing_process/MeshLocalSmoothingProcess/PerformSmoothing.md)
- [CollapseNode](../../../functions/processes/mesh_node_collapsing_process/MeshNodeCollapsingProcess/CollapseNode.md)
- [SetDistancesOnDividedElements](../../../functions/processes/parallel_distance_calculation_process/ParallelDistanceCalculationProcess<TDim>/SetDistancesOnDividedElements.md)
- [ExtendDistancesByLayer](../../../functions/processes/parallel_distance_calculation_process/ParallelDistanceCalculationProcess<TDim>/ExtendDistancesByLayer.md)
- [Execute](../../../functions/processes/tetrahedral_mesh_orientation_check/TetrahedralMeshOrientationCheck/Execute.md)
- [InternalSystemSolveWithPhysics](../../../functions/solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver/ResidualBasedBlockBuilderAndSolver/InternalSystemSolveWithPhysics.md)
- [BuildMasterSlaveConstraints](../../../functions/solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver/ResidualBasedBlockBuilderAndSolver/BuildMasterSlaveConstraints.md)
- [ApplyRHSConstraints](../../../functions/solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver_with_lagrange_multiplier/ResidualBasedBlockBuilderAndSolverWithLagrangeMultiplier/ApplyRHSConstraints.md)
- [RemoveEntities](../../../functions/sources/model_part/RemoveEntities.md)
- [ModelPartIO](../../../functions/sources/model_part_io/ModelPartIO/ModelPartIO.md)
- [ReadModelPart](../../../functions/sources/model_part_io/ModelPartIO/ReadModelPart.md)
- [WriteModelPart](../../../functions/sources/model_part_io/ModelPartIO/WriteModelPart.md)
- [ReadSubModelPartBlock](../../../functions/sources/model_part_io/ModelPartIO/ReadSubModelPartBlock.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../functions/tests/cpp_tests/processes/test_find_intersected_geometrical_objects_process/KRATOS_TEST_CASE_IN_SUITE.md)
- [testFlagSetReset](../../../functions/tests/test_flags/TestFlags/testFlagSetReset.md)
- [test_set_flag](../../../functions/tests/test_variable_utils/TestVariableUtils/test_set_flag.md)
- [RemoveElementAndBelongings](../../../functions/utilities/auxiliar_model_part_utilities/AuxiliarModelPartUtilities/RemoveElementAndBelongings.md)
- [RemoveElementsAndBelongings](../../../functions/utilities/auxiliar_model_part_utilities/AuxiliarModelPartUtilities/RemoveElementsAndBelongings.md)
- [RemoveConditionAndBelongings](../../../functions/utilities/auxiliar_model_part_utilities/AuxiliarModelPartUtilities/RemoveConditionAndBelongings.md)
- [RemoveConditionsAndBelongings](../../../functions/utilities/auxiliar_model_part_utilities/AuxiliarModelPartUtilities/RemoveConditionsAndBelongings.md)
- [EraseOuterElements](../../../functions/utilities/convect_particles_utilities/ParticleConvectUtily/EraseOuterElements.md)
- [TestIODebug](../../../functions/utilities/exact_mortar_segmentation_utility/ExactMortarIntegrationUtility<TDim, TNumNodes, TBelong, TNumNodesMaster>/TestIODebug.md)
- [ComputeTangentNodeWithLMAndSlip](../../../functions/utilities/mortar_utilities/ComputeTangentNodeWithLMAndSlip.md)
- [CalculateGradientBasedConformingExtensionBasis](../../../functions/utilities/shifted_boundary_meshless_interface_utility/ShiftedBoundaryMeshlessInterfaceUtility/CalculateGradientBasedConformingExtensionBasis.md)
- [CalculateMeshlessBasedConformingExtensionBasis](../../../functions/utilities/shifted_boundary_meshless_interface_utility/ShiftedBoundaryMeshlessInterfaceUtility/CalculateMeshlessBasedConformingExtensionBasis.md)