---
type: Python Method
title: __CheckAssignEntityDerivativesToNodesValues
resource: tests/test_sensitivity_utilities.py#L143-L153
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/tests/test_sensitivity_utilities/CalculateAnalyticalValue
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/KratosUnittest/TestCase/assertMatrixAlmostEqual
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/test_sensitivity_utilities/TestSensitivityUtilitiesTwoDimSymmetricalSquare/test_AssignElementDerivativesToNodes
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_sensitivity_utilities/TestSensitivityUtilitiesTwoDimSymmetricalSquare/test_AssignConditionDerivativesToNodes
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def __CheckAssignEntityDerivativesToNodesValues(self, nodal_neighbour_map, common_entity_map):`

# Calls

- [CalculateAnalyticalValue](../../../../functions/tests/test_sensitivity_utilities/CalculateAnalyticalValue.md)
- [assertMatrixAlmostEqual](../../../../functions/python_scripts/KratosUnittest/TestCase/assertMatrixAlmostEqual.md)

# Called by

- [test_AssignElementDerivativesToNodes](../../../../functions/tests/test_sensitivity_utilities/TestSensitivityUtilitiesTwoDimSymmetricalSquare/test_AssignElementDerivativesToNodes.md)
- [test_AssignConditionDerivativesToNodes](../../../../functions/tests/test_sensitivity_utilities/TestSensitivityUtilitiesTwoDimSymmetricalSquare/test_AssignConditionDerivativesToNodes.md)