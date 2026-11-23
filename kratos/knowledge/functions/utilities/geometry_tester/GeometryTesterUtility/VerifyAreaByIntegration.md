---
type: C++ Method
title: VerifyAreaByIntegration
resource: utilities/geometry_tester.cpp#L681-L730
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/integration/integration_point/IntegrationPoint/Weight
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/utilities/geometry_tester/GeometryTesterUtility/StreamTestTetrahedra3D4N
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/geometry_tester/GeometryTesterUtility/StreamTestTetrahedra3D10N
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/geometry_tester/GeometryTesterUtility/StreamTestTriangle2D3N
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/geometry_tester/GeometryTesterUtility/StreamTestTriangle2D6N
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/geometry_tester/GeometryTesterUtility/StreamTestQuadrilateral2D4N
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/geometry_tester/GeometryTesterUtility/StreamTestQuadrilateral2D9N
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/geometry_tester/GeometryTesterUtility/StreamTestHexahedra3D8N
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/geometry_tester/GeometryTesterUtility/StreamTestHexahedra3D20N
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/geometry_tester/GeometryTesterUtility/StreamTestHexahedra3D27N
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/geometry_tester/GeometryTesterUtility/StreamTestPrism3D6N
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/geometry_tester/GeometryTesterUtility/StreamTestPrism3D15N
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`bool GeometryTesterUtility::VerifyAreaByIntegration( GeometryType& rGeometry, GeometryType::IntegrationMethod ThisMethod, const double ReferenceArea, std::stringstream& rErrorMessage )`

# Calls

- [Weight](../../../../functions/integration/integration_point/IntegrationPoint/Weight.md)

# Called by

- [StreamTestTetrahedra3D4N](../../../../functions/utilities/geometry_tester/GeometryTesterUtility/StreamTestTetrahedra3D4N.md)
- [StreamTestTetrahedra3D10N](../../../../functions/utilities/geometry_tester/GeometryTesterUtility/StreamTestTetrahedra3D10N.md)
- [StreamTestTriangle2D3N](../../../../functions/utilities/geometry_tester/GeometryTesterUtility/StreamTestTriangle2D3N.md)
- [StreamTestTriangle2D6N](../../../../functions/utilities/geometry_tester/GeometryTesterUtility/StreamTestTriangle2D6N.md)
- [StreamTestQuadrilateral2D4N](../../../../functions/utilities/geometry_tester/GeometryTesterUtility/StreamTestQuadrilateral2D4N.md)
- [StreamTestQuadrilateral2D9N](../../../../functions/utilities/geometry_tester/GeometryTesterUtility/StreamTestQuadrilateral2D9N.md)
- [StreamTestHexahedra3D8N](../../../../functions/utilities/geometry_tester/GeometryTesterUtility/StreamTestHexahedra3D8N.md)
- [StreamTestHexahedra3D20N](../../../../functions/utilities/geometry_tester/GeometryTesterUtility/StreamTestHexahedra3D20N.md)
- [StreamTestHexahedra3D27N](../../../../functions/utilities/geometry_tester/GeometryTesterUtility/StreamTestHexahedra3D27N.md)
- [StreamTestPrism3D6N](../../../../functions/utilities/geometry_tester/GeometryTesterUtility/StreamTestPrism3D6N.md)
- [StreamTestPrism3D15N](../../../../functions/utilities/geometry_tester/GeometryTesterUtility/StreamTestPrism3D15N.md)