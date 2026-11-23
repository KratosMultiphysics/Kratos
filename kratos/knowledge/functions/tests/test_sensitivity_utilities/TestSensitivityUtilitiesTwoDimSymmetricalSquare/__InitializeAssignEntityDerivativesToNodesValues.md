---
type: Python Method
title: __InitializeAssignEntityDerivativesToNodesValues
resource: tests/test_sensitivity_utilities.py#L126-L141
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/variable_utils/SetFlag
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
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

`def __InitializeAssignEntityDerivativesToNodesValues(self, entities, nodes_per_entity):`

# Calls

- [SetFlag](../../../../functions/utilities/variable_utils/SetFlag.md)
- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)

# Called by

- [test_AssignElementDerivativesToNodes](../../../../functions/tests/test_sensitivity_utilities/TestSensitivityUtilitiesTwoDimSymmetricalSquare/test_AssignElementDerivativesToNodes.md)
- [test_AssignConditionDerivativesToNodes](../../../../functions/tests/test_sensitivity_utilities/TestSensitivityUtilitiesTwoDimSymmetricalSquare/test_AssignConditionDerivativesToNodes.md)