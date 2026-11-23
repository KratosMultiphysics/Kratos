---
type: C++ Method
title: StreamTestHexahedra3D8N
resource: utilities/geometry_tester.cpp#L390-L430
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/geometry_tester/GeometryTesterUtility/VerifyAreaByIntegration
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/geometry_tester/GeometryTesterUtility/VerifyIsInside
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/geometry_tester/GeometryTesterUtility/VerfiyShapeFunctionsValues
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/geometry_tester/GeometryTesterUtility/VerifyShapeFunctionsSecondDerivativesValues
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/utilities/geometry_tester/GeometryTesterUtility/RunTest
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/geometry_tester/TestHexahedra3D8N
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`bool GeometryTesterUtility::StreamTestHexahedra3D8N( ModelPart& rModelPart, std::stringstream& rErrorMessage )`

# Calls

- [VerifyAreaByIntegration](../../../../functions/utilities/geometry_tester/GeometryTesterUtility/VerifyAreaByIntegration.md)
- [VerifyIsInside](../../../../functions/utilities/geometry_tester/GeometryTesterUtility/VerifyIsInside.md)
- [VerfiyShapeFunctionsValues](../../../../functions/utilities/geometry_tester/GeometryTesterUtility/VerfiyShapeFunctionsValues.md)
- [VerifyShapeFunctionsSecondDerivativesValues](../../../../functions/utilities/geometry_tester/GeometryTesterUtility/VerifyShapeFunctionsSecondDerivativesValues.md)

# Called by

- [RunTest](../../../../functions/utilities/geometry_tester/GeometryTesterUtility/RunTest.md)
- [TestHexahedra3D8N](../../../../functions/utilities/geometry_tester/TestHexahedra3D8N.md)