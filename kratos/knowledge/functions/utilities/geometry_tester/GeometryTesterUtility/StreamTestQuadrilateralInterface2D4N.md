---
type: C++ Method
title: StreamTestQuadrilateralInterface2D4N
resource: utilities/geometry_tester.cpp#L361-L385
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/geometry_tester/GeometryTesterUtility/VerifyShapeFunctionsSecondDerivativesValues
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/utilities/geometry_tester/GeometryTesterUtility/RunTest
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/geometry_tester/TestQuadrilateralInterface2D4N
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`bool GeometryTesterUtility::StreamTestQuadrilateralInterface2D4N( ModelPart& rModelPart, std::stringstream& rErrorMessage )`

# Calls

- [VerifyShapeFunctionsSecondDerivativesValues](../../../../functions/utilities/geometry_tester/GeometryTesterUtility/VerifyShapeFunctionsSecondDerivativesValues.md)

# Called by

- [RunTest](../../../../functions/utilities/geometry_tester/GeometryTesterUtility/RunTest.md)
- [TestQuadrilateralInterface2D4N](../../../../functions/utilities/geometry_tester/TestQuadrilateralInterface2D4N.md)