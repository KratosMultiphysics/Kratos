---
type: Python Module
title: test_normal_utils
resource: tests/test_normal_utils.py#L1-L329
generated:
  by: okf-rs/0.7.0
relationships:
  imports:
  - target: external/os
    resolved_by: tree-sitter
    confidence: exact
  - target: external/math
    resolved_by: tree-sitter
    confidence: exact
  - target: external/kratosmultiphysics
    resolved_by: tree-sitter
    confidence: exact
  - target: external/kratosmultiphysics-kratosunittest-as-kratosunittest
    resolved_by: tree-sitter
    confidence: exact
  - target: external/kratosmultiphysics-kratos-utilities-as-kratos-utils
    resolved_by: tree-sitter
    confidence: exact
  - target: external/kratosmultiphysics-gid-output-process
    resolved_by: tree-sitter
    confidence: exact
  - target: external/kratosmultiphysics-testing-utilities
    resolved_by: tree-sitter
    confidence: exact
  member_of:
  - target: packages/KratosMultiphysics
    resolved_by: tree-sitter
    confidence: exact
---

# Contains

- [GetFilePath](../../functions/tests/test_normal_utils/GetFilePath.md)
- [PostProcess](../../functions/tests/test_normal_utils/PostProcess.md)
- [CalculateAnalyticalNormal](../../functions/tests/test_normal_utils/CalculateAnalyticalNormal.md)
- [CalculateNorm](../../functions/tests/test_normal_utils/CalculateNorm.md)
- [RemoveFiles](../../functions/tests/test_normal_utils/RemoveFiles.md)
- [FiniteDifferenceNormalShapeSensitivityTest](../../functions/tests/test_normal_utils/FiniteDifferenceNormalShapeSensitivityTest.md)
- [get_matrix_row](../../functions/tests/test_normal_utils/get_matrix_row.md)
- [check_normal_sensitivities](../../functions/tests/test_normal_utils/check_normal_sensitivities.md)
- [TestNormalUtilsCoarseSphere](../../classes/tests/test_normal_utils/TestNormalUtilsCoarseSphere.md)
- [setUpClass](../../functions/tests/test_normal_utils/TestNormalUtilsCoarseSphere/setUpClass.md)
- [tearDownClass](../../functions/tests/test_normal_utils/TestNormalUtilsCoarseSphere/tearDownClass.md)
- [setUp](../../functions/tests/test_normal_utils/TestNormalUtilsCoarseSphere/setUp.md)
- [test_ComputeSimplexNormalModelPart](../../functions/tests/test_normal_utils/TestNormalUtilsCoarseSphere/test_ComputeSimplexNormalModelPart.md)
- [test_ComputeSimplexNormalModelPartWithLineCondition](../../functions/tests/test_normal_utils/TestNormalUtilsCoarseSphere/test_ComputeSimplexNormalModelPartWithLineCondition.md)
- [test_ComputeSimplexNormalModelPartWithCustomVariable](../../functions/tests/test_normal_utils/TestNormalUtilsCoarseSphere/test_ComputeSimplexNormalModelPartWithCustomVariable.md)
- [test_ComputeSimplexNormalModelPartNonHistorical](../../functions/tests/test_normal_utils/TestNormalUtilsCoarseSphere/test_ComputeSimplexNormalModelPartNonHistorical.md)
- [test_ComputeUnitNormalModelPart](../../functions/tests/test_normal_utils/TestNormalUtilsCoarseSphere/test_ComputeUnitNormalModelPart.md)
- [test_ComputeUnitNormalModelPartWithCustomVariable](../../functions/tests/test_normal_utils/TestNormalUtilsCoarseSphere/test_ComputeUnitNormalModelPartWithCustomVariable.md)
- [test_ComputeUnitNormalModelPartNonHistorical](../../functions/tests/test_normal_utils/TestNormalUtilsCoarseSphere/test_ComputeUnitNormalModelPartNonHistorical.md)
- [test_ComputeNodesMeanNormalModelPart](../../functions/tests/test_normal_utils/TestNormalUtilsCoarseSphere/test_ComputeNodesMeanNormalModelPart.md)
- [test_ComputeNodesMeanNormalModelPartWithCustomVariable](../../functions/tests/test_normal_utils/TestNormalUtilsCoarseSphere/test_ComputeNodesMeanNormalModelPartWithCustomVariable.md)
- [test_ComputeNodesMeanNormalModelPartNonHistorical](../../functions/tests/test_normal_utils/TestNormalUtilsCoarseSphere/test_ComputeNodesMeanNormalModelPartNonHistorical.md)
- [test_InvertNormal](../../functions/tests/test_normal_utils/TestNormalUtilsCoarseSphere/test_InvertNormal.md)
- [test_CalculateShapeDerivativesOnSimplexConditions3D](../../functions/tests/test_normal_utils/TestNormalUtilsCoarseSphere/test_CalculateShapeDerivativesOnSimplexConditions3D.md)
- [TestNormalUtilsQuadSphere](../../classes/tests/test_normal_utils/TestNormalUtilsQuadSphere.md)
- [setUpClass](../../functions/tests/test_normal_utils/TestNormalUtilsQuadSphere/setUpClass.md)
- [tearDownClass](../../functions/tests/test_normal_utils/TestNormalUtilsQuadSphere/tearDownClass.md)
- [setUp](../../functions/tests/test_normal_utils/TestNormalUtilsQuadSphere/setUp.md)
- [test_ComputeUnitNormalQuadModelPart](../../functions/tests/test_normal_utils/TestNormalUtilsQuadSphere/test_ComputeUnitNormalQuadModelPart.md)
- [TestNormalUtils2DSymmetricalSquare](../../classes/tests/test_normal_utils/TestNormalUtils2DSymmetricalSquare.md)
- [setUpClass](../../functions/tests/test_normal_utils/TestNormalUtils2DSymmetricalSquare/setUpClass.md)
- [tearDownClass](../../functions/tests/test_normal_utils/TestNormalUtils2DSymmetricalSquare/tearDownClass.md)
- [setUp](../../functions/tests/test_normal_utils/TestNormalUtils2DSymmetricalSquare/setUp.md)
- [test_CalculateShapeDerivativesOnSimplexConditions2D](../../functions/tests/test_normal_utils/TestNormalUtils2DSymmetricalSquare/test_CalculateShapeDerivativesOnSimplexConditions2D.md)

# Imports

- `os`
- `math`
- `KratosMultiphysics`
- `KratosMultiphysics.KratosUnittest as KratosUnittest`
- `KratosMultiphysics.kratos_utilities as kratos_utils`
- `KratosMultiphysics.gid_output_process`
- `KratosMultiphysics.testing.utilities`

# Member of

- [KratosMultiphysics](../../packages/KratosMultiphysics.md)