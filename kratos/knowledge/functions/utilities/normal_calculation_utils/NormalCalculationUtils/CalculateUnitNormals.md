---
type: C++ Method
title: CalculateUnitNormals
resource: utilities/normal_calculation_utils.cpp#L182-L193
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/tests/test_normal_utils/TestNormalUtilsCoarseSphere/test_ComputeUnitNormalModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_normal_utils/TestNormalUtilsCoarseSphere/test_ComputeUnitNormalModelPartWithCustomVariable
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_normal_utils/TestNormalUtilsCoarseSphere/test_ComputeNodesMeanNormalModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_normal_utils/TestNormalUtilsCoarseSphere/test_ComputeNodesMeanNormalModelPartWithCustomVariable
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_normal_utils/TestNormalUtilsCoarseSphere/test_InvertNormal
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_normal_utils/TestNormalUtilsQuadSphere/test_ComputeUnitNormalQuadModelPart
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void NormalCalculationUtils::CalculateUnitNormals( ModelPart& rModelPart, const bool EnforceGenericGeometryAlgorithm, const Variable<array_1d<double,3>>& rNormalVariable )`

# Called by

- [test_ComputeUnitNormalModelPart](../../../../functions/tests/test_normal_utils/TestNormalUtilsCoarseSphere/test_ComputeUnitNormalModelPart.md)
- [test_ComputeUnitNormalModelPartWithCustomVariable](../../../../functions/tests/test_normal_utils/TestNormalUtilsCoarseSphere/test_ComputeUnitNormalModelPartWithCustomVariable.md)
- [test_ComputeNodesMeanNormalModelPart](../../../../functions/tests/test_normal_utils/TestNormalUtilsCoarseSphere/test_ComputeNodesMeanNormalModelPart.md)
- [test_ComputeNodesMeanNormalModelPartWithCustomVariable](../../../../functions/tests/test_normal_utils/TestNormalUtilsCoarseSphere/test_ComputeNodesMeanNormalModelPartWithCustomVariable.md)
- [test_InvertNormal](../../../../functions/tests/test_normal_utils/TestNormalUtilsCoarseSphere/test_InvertNormal.md)
- [test_ComputeUnitNormalQuadModelPart](../../../../functions/tests/test_normal_utils/TestNormalUtilsQuadSphere/test_ComputeUnitNormalQuadModelPart.md)