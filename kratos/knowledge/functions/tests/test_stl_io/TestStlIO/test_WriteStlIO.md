---
type: Python Method
title: test_WriteStlIO
resource: tests/test_stl_io.py#L95-L130
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_stl_io/WriteModelPartToSTL
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_stl_io/ReadModelPartFromSTL
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/communicator/Communicator/GlobalNumberOfConditions
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def test_WriteStlIO(self):`

# Calls

- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [WriteModelPartToSTL](../../../../functions/tests/test_stl_io/WriteModelPartToSTL.md)
- [ReadModelPartFromSTL](../../../../functions/tests/test_stl_io/ReadModelPartFromSTL.md)
- [GlobalNumberOfConditions](../../../../functions/sources/communicator/Communicator/GlobalNumberOfConditions.md)