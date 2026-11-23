---
type: Python Module
title: test_sensitivity_utilities
resource: tests/test_sensitivity_utilities.py#L1-L165
generated:
  by: okf-rs/0.7.0
relationships:
  imports:
  - target: external/os
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
  - target: external/kratosmultiphysics-testing-utilities
    resolved_by: tree-sitter
    confidence: exact
  member_of:
  - target: packages/KratosMultiphysics
    resolved_by: tree-sitter
    confidence: exact
---

# Contains

- [GetFilePath](../../functions/tests/test_sensitivity_utilities/GetFilePath.md)
- [RemoveFiles](../../functions/tests/test_sensitivity_utilities/RemoveFiles.md)
- [CalculateAnalyticalValue](../../functions/tests/test_sensitivity_utilities/CalculateAnalyticalValue.md)
- [CalculateEntityNeighbourMaps](../../functions/tests/test_sensitivity_utilities/CalculateEntityNeighbourMaps.md)
- [CalculateCommonEntityMap](../../functions/tests/test_sensitivity_utilities/CalculateCommonEntityMap.md)
- [TestSensitivityUtilitiesTwoDimSymmetricalSquare](../../classes/tests/test_sensitivity_utilities/TestSensitivityUtilitiesTwoDimSymmetricalSquare.md)
- [setUpClass](../../functions/tests/test_sensitivity_utilities/TestSensitivityUtilitiesTwoDimSymmetricalSquare/setUpClass.md)
- [tearDownClass](../../functions/tests/test_sensitivity_utilities/TestSensitivityUtilitiesTwoDimSymmetricalSquare/tearDownClass.md)
- [setUp](../../functions/tests/test_sensitivity_utilities/TestSensitivityUtilitiesTwoDimSymmetricalSquare/setUp.md)
- [test_AssignElementDerivativesToNodes](../../functions/tests/test_sensitivity_utilities/TestSensitivityUtilitiesTwoDimSymmetricalSquare/test_AssignElementDerivativesToNodes.md)
- [test_AssignConditionDerivativesToNodes](../../functions/tests/test_sensitivity_utilities/TestSensitivityUtilitiesTwoDimSymmetricalSquare/test_AssignConditionDerivativesToNodes.md)
- [__InitializeAssignEntityDerivativesToNodesValues](../../functions/tests/test_sensitivity_utilities/TestSensitivityUtilitiesTwoDimSymmetricalSquare/__InitializeAssignEntityDerivativesToNodesValues.md)
- [__CheckAssignEntityDerivativesToNodesValues](../../functions/tests/test_sensitivity_utilities/TestSensitivityUtilitiesTwoDimSymmetricalSquare/__CheckAssignEntityDerivativesToNodesValues.md)
- [TestSensitivityUtilitiesGetSensitivityVariableName](../../classes/tests/test_sensitivity_utilities/TestSensitivityUtilitiesGetSensitivityVariableName.md)
- [test_GetSensitivityVariableName](../../functions/tests/test_sensitivity_utilities/TestSensitivityUtilitiesGetSensitivityVariableName/test_GetSensitivityVariableName.md)

# Imports

- `os`
- `KratosMultiphysics`
- `KratosMultiphysics.KratosUnittest as KratosUnittest`
- `KratosMultiphysics.kratos_utilities as kratos_utils`
- `KratosMultiphysics.testing.utilities`

# Member of

- [KratosMultiphysics](../../packages/KratosMultiphysics.md)