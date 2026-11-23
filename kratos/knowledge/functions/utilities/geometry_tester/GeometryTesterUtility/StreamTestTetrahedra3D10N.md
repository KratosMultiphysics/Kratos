---
type: C++ Method
title: StreamTestTetrahedra3D10N
resource: utilities/geometry_tester.cpp#L142-L179
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/geometry_tester/GeometryTesterUtility/VerifyAreaByIntegration
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/utilities/geometry_tester/GeometryTesterUtility/RunTest
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/geometry_tester/TestTetrahedra3D10N
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`bool GeometryTesterUtility::StreamTestTetrahedra3D10N( ModelPart& rModelPart, std::stringstream& rErrorMessage )`

# Calls

- [VerifyAreaByIntegration](../../../../functions/utilities/geometry_tester/GeometryTesterUtility/VerifyAreaByIntegration.md)

# Called by

- [RunTest](../../../../functions/utilities/geometry_tester/GeometryTesterUtility/RunTest.md)
- [TestTetrahedra3D10N](../../../../functions/utilities/geometry_tester/TestTetrahedra3D10N.md)