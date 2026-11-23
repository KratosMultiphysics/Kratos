---
type: C++ Method
title: StreamTestTriangle2D6N
resource: utilities/geometry_tester.cpp#L231-L268
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/geometry_tester/GeometryTesterUtility/VerifyAreaByIntegration
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/geometry_tester/GeometryTesterUtility/VerifyShapeFunctionsSecondDerivativesValues
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/utilities/geometry_tester/GeometryTesterUtility/RunTest
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/geometry_tester/TestTriangle2D6N
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`bool GeometryTesterUtility::StreamTestTriangle2D6N( ModelPart& rModelPart, std::stringstream& rErrorMessage )`

# Calls

- [VerifyAreaByIntegration](../../../../functions/utilities/geometry_tester/GeometryTesterUtility/VerifyAreaByIntegration.md)
- [VerifyShapeFunctionsSecondDerivativesValues](../../../../functions/utilities/geometry_tester/GeometryTesterUtility/VerifyShapeFunctionsSecondDerivativesValues.md)

# Called by

- [RunTest](../../../../functions/utilities/geometry_tester/GeometryTesterUtility/RunTest.md)
- [TestTriangle2D6N](../../../../functions/utilities/geometry_tester/TestTriangle2D6N.md)