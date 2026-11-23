---
type: Python Method
title: test_TestVoxelMeshModelerQuadraticElementInput
resource: tests/test_voxel_modeler.py#L143-L199
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/GetString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/KratosUnittest/TestCase/assertEqualTolerance
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def test_TestVoxelMeshModelerQuadraticElementInput(self): ## create model`

# Calls

- [GetString](../../../../functions/sources/kratos_parameters/Parameters/GetString.md)
- [assertEqualTolerance](../../../../functions/python_scripts/KratosUnittest/TestCase/assertEqualTolerance.md)
- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)