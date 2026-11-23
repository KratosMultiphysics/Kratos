---
type: C++ Function
title: SetHistoricalVariableToZero
resource: utilities/variable_utils.h#L757-L764
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/processes/calculate_distance_to_path_process/CalculateDistanceToPathProcess<THistorical>/Execute
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/ResidualBasedAdjointBossakScheme/SetToZero_AdjointVars
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_coordinate_transformation_utils/TestCoordinateTransformationUtilitiesCoarseSphere/setUp
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_coordinate_transformation_utils/TestCoordinateTransformationUtilities2DSymmetricalSquare/setUp
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_normal_utils/TestNormalUtilsCoarseSphere/setUp
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_normal_utils/TestNormalUtilsQuadSphere/setUp
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_normal_utils/TestNormalUtils2DSymmetricalSquare/setUp
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_variable_utils/TestVariableUtils/test_set_nodal_historical_variable_to_zero
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void SetHistoricalVariableToZero( const Variable< TType >& rVariable, NodesContainerType& rNodes)`

# Called by

- [Execute](../../../functions/processes/calculate_distance_to_path_process/CalculateDistanceToPathProcess<THistorical>/Execute.md)
- [SetToZero_AdjointVars](../../../functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/ResidualBasedAdjointBossakScheme/SetToZero_AdjointVars.md)
- [setUp](../../../functions/tests/test_coordinate_transformation_utils/TestCoordinateTransformationUtilitiesCoarseSphere/setUp.md)
- [setUp](../../../functions/tests/test_coordinate_transformation_utils/TestCoordinateTransformationUtilities2DSymmetricalSquare/setUp.md)
- [setUp](../../../functions/tests/test_normal_utils/TestNormalUtilsCoarseSphere/setUp.md)
- [setUp](../../../functions/tests/test_normal_utils/TestNormalUtilsQuadSphere/setUp.md)
- [setUp](../../../functions/tests/test_normal_utils/TestNormalUtils2DSymmetricalSquare/setUp.md)
- [test_set_nodal_historical_variable_to_zero](../../../functions/tests/test_variable_utils/TestVariableUtils/test_set_nodal_historical_variable_to_zero.md)