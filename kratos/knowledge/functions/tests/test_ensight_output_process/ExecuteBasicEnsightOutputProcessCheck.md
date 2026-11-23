---
type: Python Function
title: ExecuteBasicEnsightOutputProcessCheck
resource: tests/test_ensight_output_process.py#L197-L299
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/SetString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_ensight_output_process/SetupEnsightOutputProcess
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetStringArray
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/test_ensight_output_process/TestEnsightOutputProcess/test_ascii_ensight_output_2D
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_ensight_output_process/TestEnsightOutputProcess/test_ascii_ensight_output_3D
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_ensight_output_process/TestEnsightOutputProcess/test_ascii_ensight_output_quad_3D
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_ensight_output_process/TestEnsightOutputProcess/test_ascii_ensight_output_quad_prism_3D
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_ensight_output_process/TestEnsightOutputProcess/test_ascii_ensight_output_quad_hexahedra_3D
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_ensight_output_process/TestEnsightOutputProcess/test_binary_ensight_output_2D
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_ensight_output_process/TestEnsightOutputProcess/test_binary_ensight_output_3D
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_ensight_output_process/TestEnsightOutputProcess/test_binary_ensight_output_quad_3D
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_ensight_output_process/TestEnsightOutputProcess/test_binary_ensight_output_quad_prism_3D
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_ensight_output_process/TestEnsightOutputProcess/test_binary_ensight_output_quad_hexahedra_3D
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def ExecuteBasicEnsightOutputProcessCheck(file_format = "ascii", setup = "2D"):`

# Calls

- [SetString](../../../functions/sources/kratos_parameters/Parameters/SetString.md)
- [SetupEnsightOutputProcess](../../../functions/tests/test_ensight_output_process/SetupEnsightOutputProcess.md)
- [GetStringArray](../../../functions/sources/kratos_parameters/Parameters/GetStringArray.md)

# Called by

- [test_ascii_ensight_output_2D](../../../functions/tests/test_ensight_output_process/TestEnsightOutputProcess/test_ascii_ensight_output_2D.md)
- [test_ascii_ensight_output_3D](../../../functions/tests/test_ensight_output_process/TestEnsightOutputProcess/test_ascii_ensight_output_3D.md)
- [test_ascii_ensight_output_quad_3D](../../../functions/tests/test_ensight_output_process/TestEnsightOutputProcess/test_ascii_ensight_output_quad_3D.md)
- [test_ascii_ensight_output_quad_prism_3D](../../../functions/tests/test_ensight_output_process/TestEnsightOutputProcess/test_ascii_ensight_output_quad_prism_3D.md)
- [test_ascii_ensight_output_quad_hexahedra_3D](../../../functions/tests/test_ensight_output_process/TestEnsightOutputProcess/test_ascii_ensight_output_quad_hexahedra_3D.md)
- [test_binary_ensight_output_2D](../../../functions/tests/test_ensight_output_process/TestEnsightOutputProcess/test_binary_ensight_output_2D.md)
- [test_binary_ensight_output_3D](../../../functions/tests/test_ensight_output_process/TestEnsightOutputProcess/test_binary_ensight_output_3D.md)
- [test_binary_ensight_output_quad_3D](../../../functions/tests/test_ensight_output_process/TestEnsightOutputProcess/test_binary_ensight_output_quad_3D.md)
- [test_binary_ensight_output_quad_prism_3D](../../../functions/tests/test_ensight_output_process/TestEnsightOutputProcess/test_binary_ensight_output_quad_prism_3D.md)
- [test_binary_ensight_output_quad_hexahedra_3D](../../../functions/tests/test_ensight_output_process/TestEnsightOutputProcess/test_binary_ensight_output_quad_hexahedra_3D.md)