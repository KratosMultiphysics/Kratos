---
type: C++ Function
title: VectorizedFindHelper
resource: python/add_geometrical_utilities_to_python.cpp#L48-L87
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/parallel_utilities/IndexPartition/IndexPartition
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/binbased_fast_point_locator/FindPointOnMeshSimplified
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`std::tuple<pybind11::array_t<int>, pybind11::array_t<int>, pybind11::array_t<double> > VectorizedFindHelper(BinBasedFastPointLocator < TDim >& rSelf, const Matrix& rCoords )`

# Calls

- [IndexPartition](../../../functions/utilities/parallel_utilities/IndexPartition/IndexPartition.md)
- [FindPointOnMeshSimplified](../../../functions/utilities/binbased_fast_point_locator/FindPointOnMeshSimplified.md)
- [GetGeometry](../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)