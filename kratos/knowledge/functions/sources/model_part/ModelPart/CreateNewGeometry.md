---
type: C++ Method
title: CreateNewGeometry
resource: sources/model_part.cpp#L1618-L1635
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/model_part/IsSubModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/geometry/PointsNumber
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`ModelPart::GeometryType::Pointer ModelPart::CreateNewGeometry( const std::string& rGeometryTypeName, const std::vector<IndexType>& rGeometryNodeIds )`

# Calls

- [IsSubModelPart](../../../../functions/includes/model_part/IsSubModelPart.md)
- [PointsNumber](../../../../functions/geometries/geometry/PointsNumber.md)