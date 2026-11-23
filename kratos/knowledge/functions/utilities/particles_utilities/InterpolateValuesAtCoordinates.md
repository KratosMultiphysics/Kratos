---
type: C++ Function
title: InterpolateValuesAtCoordinates
resource: utilities/particles_utilities.h#L243-L280
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/parallel_utilities/IndexPartition/IndexPartition
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static std::pair< DenseVector<bool>, std::vector<TDataType> > InterpolateValuesAtCoordinates( BinBasedFastPointLocator<TDim>& rLocator, const Matrix& rCoordinates, const Variable<TDataType>& rInterpolationVariable, const double SearchTolerance )`

# Calls

- [IndexPartition](../../../functions/utilities/parallel_utilities/IndexPartition/IndexPartition.md)
- [GetGeometry](../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)