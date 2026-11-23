---
type: C++ Method
title: AssignEntityDerivativesToNodes
resource: utilities/sensitivity_utilities.cpp#L171-L313
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/variable_utils/SetFlag
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/geometry/PointsNumber
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/sensitivity_utilities/SensitivityUtilities/ComputeEntityGeometryNeighbourNodeMap
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/sensitivity_utilities/SensitivityUtilities/GetMatrixSubBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/sensitivity_utilities/SensitivityUtilities/AddMatrixSubBlock
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void SensitivityUtilities::AssignEntityDerivativesToNodes( ModelPart& rModelPart, const int DerivativeDimension, const Variable<Matrix>& rDerivativeVariable, const std::unordered_map<int, std::vector<int>>& rNeighbourNodeIdsMap, const double Weight, const Flags& rFlag, const bool CheckValue)`

# Calls

- [SetFlag](../../../../functions/utilities/variable_utils/SetFlag.md)
- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [PointsNumber](../../../../functions/geometries/geometry/PointsNumber.md)
- [ComputeEntityGeometryNeighbourNodeMap](../../../../functions/utilities/sensitivity_utilities/SensitivityUtilities/ComputeEntityGeometryNeighbourNodeMap.md)
- [GetMatrixSubBlock](../../../../functions/utilities/sensitivity_utilities/SensitivityUtilities/GetMatrixSubBlock.md)
- [AddMatrixSubBlock](../../../../functions/utilities/sensitivity_utilities/SensitivityUtilities/AddMatrixSubBlock.md)