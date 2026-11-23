---
type: C++ Method
title: ComputeElementSkin
resource: utilities/embedded_skin_utility.cpp#L118-L177
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/embedded_skin_utility/EmbeddedSkinUtility<TDim>/SetDivideGeometryUtility
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/divide_geometry/DivideGeometry<TPointType>/GetPositiveInterfaces
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/geometry/PointsNumber
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/embedded_skin_utility/EmbeddedSkinUtility<TDim>/pCreateNewCondition
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/utilities/embedded_skin_utility/EmbeddedSkinUtility<TDim>/GenerateSkin
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void EmbeddedSkinUtility<TDim>::ComputeElementSkin( const Element::Pointer pElement, const Vector &rNodalDistances, unsigned int &rTempNodeId, unsigned int &rTempCondId, Properties::Pointer pCondProp, std::vector<Node::Pointer> &rNewNodesVect, std::vector<Condition::Pointer> &rNewCondsVect)`

# Calls

- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [SetDivideGeometryUtility](../../../../functions/utilities/embedded_skin_utility/EmbeddedSkinUtility<TDim>/SetDivideGeometryUtility.md)
- [GetPositiveInterfaces](../../../../functions/utilities/divide_geometry/DivideGeometry<TPointType>/GetPositiveInterfaces.md)
- [PointsNumber](../../../../functions/geometries/geometry/PointsNumber.md)
- [pCreateNewCondition](../../../../functions/utilities/embedded_skin_utility/EmbeddedSkinUtility<TDim>/pCreateNewCondition.md)

# Called by

- [GenerateSkin](../../../../functions/utilities/embedded_skin_utility/EmbeddedSkinUtility<TDim>/GenerateSkin.md)