---
type: C++ Method
title: ComputeEntityGeometryNeighbourNodeMap
resource: utilities/sensitivity_utilities.cpp#L104-L154
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/geometries/geometry/PointsNumber
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/utilities/sensitivity_utilities/SensitivityUtilities/AssignEntityDerivativesToNodes
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void SensitivityUtilities::ComputeEntityGeometryNeighbourNodeMap( std::unordered_map<int, std::unordered_map<int, int>>& rDerivativeNodesMap, const std::unordered_map<int, std::vector<int>>& rNeighbourNodeIdsMap, const Geometry<ModelPart::NodeType>& rEntityGeometry, const Flags& rFlag, const bool CheckValue)`

# Calls

- [PointsNumber](../../../../functions/geometries/geometry/PointsNumber.md)

# Called by

- [AssignEntityDerivativesToNodes](../../../../functions/utilities/sensitivity_utilities/SensitivityUtilities/AssignEntityDerivativesToNodes.md)