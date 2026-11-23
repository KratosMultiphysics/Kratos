---
type: C++ Method
title: GlobalSpaceDerivatives
resource: geometries/nurbs_curve_on_surface_geometry.h#L863-L906
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/geometries/nurbs_shape_function_utilities/nurbs_utilities/GetBinomCoefficient
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void GlobalSpaceDerivatives( std::vector<CoordinatesArrayType>& rGlobalSpaceDerivatives, const CoordinatesArrayType& rCoordinates, const SizeType DerivativeOrder) const override`

# Calls

- [GetBinomCoefficient](../../../../functions/geometries/nurbs_shape_function_utilities/nurbs_utilities/GetBinomCoefficient.md)