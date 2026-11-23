---
type: Python Method
title: test_WriteObjIO
resource: tests/test_obj_io.py#L242-L294
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/variable_utils/SetNonHistoricalVariableToZero
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_obj_io/WriteModelPartToOBJ
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_obj_io/ReadModelPartFromOBJ
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_obj_io/calculate_analytical_normal
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_obj_io/calculate_norm
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def test_WriteObjIO(self):`

# Calls

- [SetNonHistoricalVariableToZero](../../../../functions/utilities/variable_utils/SetNonHistoricalVariableToZero.md)
- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [WriteModelPartToOBJ](../../../../functions/tests/test_obj_io/WriteModelPartToOBJ.md)
- [ReadModelPartFromOBJ](../../../../functions/tests/test_obj_io/ReadModelPartFromOBJ.md)
- [calculate_analytical_normal](../../../../functions/tests/test_obj_io/calculate_analytical_normal.md)
- [calculate_norm](../../../../functions/tests/test_obj_io/calculate_norm.md)