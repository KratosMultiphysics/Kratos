---
type: C++ Function
title: InterpolateMeshVariableToSkinSpecialization
resource: utilities/embedded_skin_utility.h#L271-L322
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/node/SolutionStepsDataHas
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/embedded_skin_utility/EmbeddedSkinUtility<TDim>/SetDistancesVector
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/embedded_skin_utility/EmbeddedSkinUtility<TDim>/pCreateModifiedShapeFunctions
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/embedded_skin_utility/EmbeddedSkinUtility<TDim>/GetModifiedShapeFunctionsValuesOnEdge
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/geometry/PointsNumber
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void InterpolateMeshVariableToSkinSpecialization( const Variable<TVarType> &rMeshVariable, const Variable<TVarType> &rSkinVariable, const std::string &rInterfaceSide = "positive")`

# Calls

- [SolutionStepsDataHas](../../../functions/includes/node/SolutionStepsDataHas.md)
- [SetDistancesVector](../../../functions/utilities/embedded_skin_utility/EmbeddedSkinUtility<TDim>/SetDistancesVector.md)
- [pCreateModifiedShapeFunctions](../../../functions/utilities/embedded_skin_utility/EmbeddedSkinUtility<TDim>/pCreateModifiedShapeFunctions.md)
- [GetModifiedShapeFunctionsValuesOnEdge](../../../functions/utilities/embedded_skin_utility/EmbeddedSkinUtility<TDim>/GetModifiedShapeFunctionsValuesOnEdge.md)
- [GetGeometry](../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [PointsNumber](../../../functions/geometries/geometry/PointsNumber.md)