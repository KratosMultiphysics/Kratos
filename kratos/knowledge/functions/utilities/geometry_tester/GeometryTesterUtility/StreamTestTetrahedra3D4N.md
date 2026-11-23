---
type: C++ Method
title: StreamTestTetrahedra3D4N
resource: utilities/geometry_tester.cpp#L96-L137
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
  - target: functions/utilities/geometry_tester/TestTetrahedra3D4N
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`bool GeometryTesterUtility::StreamTestTetrahedra3D4N( ModelPart& rModelPart, std::stringstream& rErrorMessage )`

# Calls

- [VerifyAreaByIntegration](../../../../functions/utilities/geometry_tester/GeometryTesterUtility/VerifyAreaByIntegration.md)
- [VerifyShapeFunctionsSecondDerivativesValues](../../../../functions/utilities/geometry_tester/GeometryTesterUtility/VerifyShapeFunctionsSecondDerivativesValues.md)

# Called by

- [RunTest](../../../../functions/utilities/geometry_tester/GeometryTesterUtility/RunTest.md)
- [TestTetrahedra3D4N](../../../../functions/utilities/geometry_tester/TestTetrahedra3D4N.md)