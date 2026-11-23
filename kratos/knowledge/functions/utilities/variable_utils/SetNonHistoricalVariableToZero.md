---
type: C++ Function
title: SetNonHistoricalVariableToZero
resource: utilities/variable_utils.h#L724-L731
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/processes/calculate_discontinuous_distance_to_skin_process/CalculateEmbeddedVariableFromSkinSpecialization
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_distance_to_path_process/CalculateDistanceToPathProcess<THistorical>/Execute
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_nodal_distance_to_skin_process/CalculateNodalDistanceToSkinProcess/Execute
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/edge_based_gradient_recovery_process/EdgeBasedGradientRecoveryProcess/ExecuteInitialize
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/simple_mortar_mapper_process/SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>/ResetNodalArea
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/time_averaging_process/TimeAveragingProcess/ExecuteInitialize
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/ResidualBasedAdjointBossakScheme/Initialize
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/ResidualBasedAdjointBossakScheme/CalculateNodeNeighbourCount
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_calculate_nodal_distance_to_skin_process/TestCalculateNodalDistanceToSkinProcessCoarseSphere/test_ComputeDistanceToSkin
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_obj_io/TestObjIO/test_WriteObjIO
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_variable_utils/TestVariableUtils/test_set_variable_to_zero
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void SetNonHistoricalVariableToZero( const Variable< TType >& rVariable, TContainerType& rContainer)`

# Called by

- [CalculateEmbeddedVariableFromSkinSpecialization](../../../functions/processes/calculate_discontinuous_distance_to_skin_process/CalculateEmbeddedVariableFromSkinSpecialization.md)
- [Execute](../../../functions/processes/calculate_distance_to_path_process/CalculateDistanceToPathProcess<THistorical>/Execute.md)
- [Execute](../../../functions/processes/calculate_nodal_distance_to_skin_process/CalculateNodalDistanceToSkinProcess/Execute.md)
- [ExecuteInitialize](../../../functions/processes/edge_based_gradient_recovery_process/EdgeBasedGradientRecoveryProcess/ExecuteInitialize.md)
- [ResetNodalArea](../../../functions/processes/simple_mortar_mapper_process/SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>/ResetNodalArea.md)
- [ExecuteInitialize](../../../functions/processes/time_averaging_process/TimeAveragingProcess/ExecuteInitialize.md)
- [Initialize](../../../functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/ResidualBasedAdjointBossakScheme/Initialize.md)
- [CalculateNodeNeighbourCount](../../../functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/ResidualBasedAdjointBossakScheme/CalculateNodeNeighbourCount.md)
- [test_ComputeDistanceToSkin](../../../functions/tests/test_calculate_nodal_distance_to_skin_process/TestCalculateNodalDistanceToSkinProcessCoarseSphere/test_ComputeDistanceToSkin.md)
- [test_WriteObjIO](../../../functions/tests/test_obj_io/TestObjIO/test_WriteObjIO.md)
- [test_set_variable_to_zero](../../../functions/tests/test_variable_utils/TestVariableUtils/test_set_variable_to_zero.md)