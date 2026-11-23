---
type: C++ Method
title: RecursivelyValidateAndAssignDefaults
resource: sources/kratos_parameters.cpp#L1428-L1436
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/RecursivelyAddMissingParameters
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/processes/apply_periodic_boundary_condition_process/ApplyPeriodicConditionProcess/ApplyPeriodicConditionProcess
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_discontinuous_distance_to_skin_process/CalculateDiscontinuousDistanceToSkinProcess<TDim>/CalculateDiscontinuousDistanceToSkinProcess
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_nodal_area_process/CalculateNodalAreaProcess<THistorical>/CalculateNodalAreaProcess
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/check_same_modelpart_using_skin_distance_process/CheckSameModelPartUsingSkinDistanceProcess
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/compute_nodal_gradient_process/ComputeNodalGradientProcess<TOutputHistorical>/PrepareMemberVariablesFromSettings
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/entity_erase_process/EntitiesEraseProcess<TEntity>/EntitiesEraseProcess
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/find_intersected_geometrical_objects_process/FindIntersectedGeometricalObjectsProcess/FindIntersectedGeometricalObjectsProcess
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/find_intersected_geometrical_objects_with_obb_process/FindIntersectedGeometricalObjectsWithOBBProcess/FindIntersectedGeometricalObjectsWithOBBProcess
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/simple_mortar_mapper_process/SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>/SimpleMortarMapperProcess
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver_with_lagrange_multiplier/ResidualBasedBlockBuilderAndSolverWithLagrangeMultiplier/ValidateAndAssignParameters
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/spatial_containers/specialized_spatial_search/SpecializedSpatialSearch
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/spatial_containers/specialized_spatial_search_factory/SpecializedSpatialSearchFactory/SpecializedSpatialSearchFactory
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/sources/test_kratos_parameters/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_kratos_parameters/TestParameters/test_recursive_validation_fails_error_on_first_level
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_kratos_parameters/TestParameters/test_recursive_validation_4_levels
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/read_materials_utility/ReadMaterialsUtility/ReadMaterialsUtility
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/svd_utils/SingularValueDecomposition
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void Parameters::RecursivelyValidateAndAssignDefaults(const Parameters& rDefaultParameters)`

# Calls

- [RecursivelyAddMissingParameters](../../../../functions/sources/kratos_parameters/Parameters/RecursivelyAddMissingParameters.md)

# Called by

- [ApplyPeriodicConditionProcess](../../../../functions/processes/apply_periodic_boundary_condition_process/ApplyPeriodicConditionProcess/ApplyPeriodicConditionProcess.md)
- [CalculateDiscontinuousDistanceToSkinProcess](../../../../functions/processes/calculate_discontinuous_distance_to_skin_process/CalculateDiscontinuousDistanceToSkinProcess<TDim>/CalculateDiscontinuousDistanceToSkinProcess.md)
- [CalculateNodalAreaProcess](../../../../functions/processes/calculate_nodal_area_process/CalculateNodalAreaProcess<THistorical>/CalculateNodalAreaProcess.md)
- [CheckSameModelPartUsingSkinDistanceProcess](../../../../functions/processes/check_same_modelpart_using_skin_distance_process/CheckSameModelPartUsingSkinDistanceProcess.md)
- [PrepareMemberVariablesFromSettings](../../../../functions/processes/compute_nodal_gradient_process/ComputeNodalGradientProcess<TOutputHistorical>/PrepareMemberVariablesFromSettings.md)
- [EntitiesEraseProcess](../../../../functions/processes/entity_erase_process/EntitiesEraseProcess<TEntity>/EntitiesEraseProcess.md)
- [FindIntersectedGeometricalObjectsProcess](../../../../functions/processes/find_intersected_geometrical_objects_process/FindIntersectedGeometricalObjectsProcess/FindIntersectedGeometricalObjectsProcess.md)
- [FindIntersectedGeometricalObjectsWithOBBProcess](../../../../functions/processes/find_intersected_geometrical_objects_with_obb_process/FindIntersectedGeometricalObjectsWithOBBProcess/FindIntersectedGeometricalObjectsWithOBBProcess.md)
- [SimpleMortarMapperProcess](../../../../functions/processes/simple_mortar_mapper_process/SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>/SimpleMortarMapperProcess.md)
- [ValidateAndAssignParameters](../../../../functions/solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver_with_lagrange_multiplier/ResidualBasedBlockBuilderAndSolverWithLagrangeMultiplier/ValidateAndAssignParameters.md)
- [SpecializedSpatialSearch](../../../../functions/spatial_containers/specialized_spatial_search/SpecializedSpatialSearch.md)
- [SpecializedSpatialSearchFactory](../../../../functions/spatial_containers/specialized_spatial_search_factory/SpecializedSpatialSearchFactory/SpecializedSpatialSearchFactory.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/sources/test_kratos_parameters/KRATOS_TEST_CASE_IN_SUITE.md)
- [test_recursive_validation_fails_error_on_first_level](../../../../functions/tests/test_kratos_parameters/TestParameters/test_recursive_validation_fails_error_on_first_level.md)
- [test_recursive_validation_4_levels](../../../../functions/tests/test_kratos_parameters/TestParameters/test_recursive_validation_4_levels.md)
- [ReadMaterialsUtility](../../../../functions/utilities/read_materials_utility/ReadMaterialsUtility/ReadMaterialsUtility.md)
- [SingularValueDecomposition](../../../../functions/utilities/svd_utils/SingularValueDecomposition.md)