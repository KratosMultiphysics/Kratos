---
type: C++ Method
title: VerfiyShapeFunctionsValues
resource: utilities/geometry_tester.cpp#L934-L959
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/amatrix_interface/norm_2
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/utilities/geometry_tester/GeometryTesterUtility/StreamTestHexahedra3D8N
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`bool GeometryTesterUtility::VerfiyShapeFunctionsValues( GeometryType& rGeometry, GeometryType::CoordinatesArrayType& rGlobalCoordinates, std::stringstream& rErrorMessage )`

# Calls

- [norm_2](../../../../functions/includes/amatrix_interface/norm_2.md)

# Called by

- [StreamTestHexahedra3D8N](../../../../functions/utilities/geometry_tester/GeometryTesterUtility/StreamTestHexahedra3D8N.md)