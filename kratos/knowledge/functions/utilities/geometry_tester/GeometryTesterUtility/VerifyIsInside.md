---
type: C++ Method
title: VerifyIsInside
resource: utilities/geometry_tester.cpp#L914-L929
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/utilities/geometry_tester/GeometryTesterUtility/StreamTestHexahedra3D8N
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`bool GeometryTesterUtility::VerifyIsInside( GeometryType& rGeometry, GeometryType::CoordinatesArrayType& rGlobalCoordinates, bool ExpectedResult, std::stringstream& rErrorMessage)`

# Called by

- [StreamTestHexahedra3D8N](../../../../functions/utilities/geometry_tester/GeometryTesterUtility/StreamTestHexahedra3D8N.md)