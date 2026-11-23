---
type: C++ Function
title: ReadNurbsCurve
resource: input_output/cad_json_input.h#L635-L687
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/GetBool
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetInt
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/cad_json_input/CadJsonInput/ReadControlPointWeightVector
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static typename NurbsCurveGeometry<TWorkingSpaceDimension, PointerVector<TThisNodeType>>::Pointer ReadNurbsCurve( const Parameters rParameters, ModelPart& rModelPart, ProjectionAlgorithm ProjectionAlgorithmType = ProjectionAlgorithm::NewtonRaphson, SizeType EchoLevel = 0)`

# Calls

- [GetBool](../../../functions/sources/kratos_parameters/Parameters/GetBool.md)
- [GetInt](../../../functions/sources/kratos_parameters/Parameters/GetInt.md)
- [ReadControlPointWeightVector](../../../functions/input_output/cad_json_input/CadJsonInput/ReadControlPointWeightVector.md)