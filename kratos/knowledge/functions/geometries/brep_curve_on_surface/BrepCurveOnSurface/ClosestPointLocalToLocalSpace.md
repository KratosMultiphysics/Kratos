---
type: C++ Method
title: ClosestPointLocalToLocalSpace
resource: geometries/brep_curve_on_surface.h#L386-L412
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/geometries/nurbs_shape_function_utilities/nurbs_interval/NurbsInterval/MinParameter
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/nurbs_shape_function_utilities/nurbs_interval/NurbsInterval/MaxParameter
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`virtual int ClosestPointLocalToLocalSpace( const CoordinatesArrayType& rPointLocalCoordinates, CoordinatesArrayType& rClosestPointLocalCoordinates, const double Tolerance = std::numeric_limits<double>::epsilon() ) const override`

# Calls

- [MinParameter](../../../../functions/geometries/nurbs_shape_function_utilities/nurbs_interval/NurbsInterval/MinParameter.md)
- [MaxParameter](../../../../functions/geometries/nurbs_shape_function_utilities/nurbs_interval/NurbsInterval/MaxParameter.md)