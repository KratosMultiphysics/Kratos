---
type: C++ Function
title: FindPointOnMesh
resource: utilities/binbased_fast_point_locator.h#L184-L218
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`bool FindPointOnMesh( const array_1d<double, 3 >& rCoordinates, array_1d<double, TDim + 1 >& rNShapeFunction, typename EntityType::Pointer& pEntity, ResultIteratorType ItResultBegin, const SizeType MaxNumberOfResults = 1000, const double Tolerance = 1.0e-5 )`

# Calls

- [GetGeometry](../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)