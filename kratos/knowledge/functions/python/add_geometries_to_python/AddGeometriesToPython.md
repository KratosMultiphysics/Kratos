---
type: C++ Function
title: AddGeometriesToPython
resource: python/add_geometries_to_python.cpp#L116-L508
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/geometries/geometry/GetDefaultIntegrationMethod
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/geometry/PointsNumber
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/brep_surface/BrepSurface/KnotsU
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/brep_surface/BrepSurface/KnotsV
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/brep_surface/BrepSurface/ComputeSpanTriangulationLocalSpace
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/python/kratos_python/PYBIND11_MODULE
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void AddGeometriesToPython(pybind11::module& m)`

# Calls

- [GetDefaultIntegrationMethod](../../../functions/geometries/geometry/GetDefaultIntegrationMethod.md)
- [PointsNumber](../../../functions/geometries/geometry/PointsNumber.md)
- [KnotsU](../../../functions/geometries/brep_surface/BrepSurface/KnotsU.md)
- [KnotsV](../../../functions/geometries/brep_surface/BrepSurface/KnotsV.md)
- [ComputeSpanTriangulationLocalSpace](../../../functions/geometries/brep_surface/BrepSurface/ComputeSpanTriangulationLocalSpace.md)

# Called by

- [PYBIND11_MODULE](../../../functions/python/kratos_python/PYBIND11_MODULE.md)