---
type: Python Function
title: ExecuteEnsightTensorOutputProcessCheck
resource: tests/test_ensight_output_process.py#L348-L421
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/tests/test_ensight_output_process/SetupEnsightOutputProcess
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_ensight_output_process/SetSolutionWithMatrixVariables
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_ensight_output_process/CheckFileNonEmpty
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/test_ensight_output_process/TestEnsightOutputProcess/test_ascii_ensight6_symmetric_tensor_2D
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_ensight_output_process/TestEnsightOutputProcess/test_binary_ensight6_symmetric_tensor_2D
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_ensight_output_process/TestEnsightOutputProcess/test_ascii_ensightgold_symmetric_tensor_2D
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_ensight_output_process/TestEnsightOutputProcess/test_binary_ensightgold_symmetric_tensor_2D
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_ensight_output_process/TestEnsightOutputProcess/test_ascii_ensight6_symmetric_tensor_3D
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_ensight_output_process/TestEnsightOutputProcess/test_binary_ensight6_symmetric_tensor_3D
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_ensight_output_process/TestEnsightOutputProcess/test_ascii_ensightgold_symmetric_tensor_3D
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_ensight_output_process/TestEnsightOutputProcess/test_binary_ensightgold_symmetric_tensor_3D
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_ensight_output_process/TestEnsightOutputProcess/test_ascii_ensight6_asymmetric_tensor_2D
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_ensight_output_process/TestEnsightOutputProcess/test_binary_ensight6_asymmetric_tensor_2D
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_ensight_output_process/TestEnsightOutputProcess/test_ascii_ensightgold_asymmetric_tensor_2D
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_ensight_output_process/TestEnsightOutputProcess/test_binary_ensightgold_asymmetric_tensor_2D
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_ensight_output_process/TestEnsightOutputProcess/test_ascii_ensight6_asymmetric_tensor_3D
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_ensight_output_process/TestEnsightOutputProcess/test_binary_ensight6_asymmetric_tensor_3D
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_ensight_output_process/TestEnsightOutputProcess/test_ascii_ensightgold_asymmetric_tensor_3D
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_ensight_output_process/TestEnsightOutputProcess/test_binary_ensightgold_asymmetric_tensor_3D
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def ExecuteEnsightTensorOutputProcessCheck(file_format="ascii", ensight_format="gold", setup="2D", tensor_type="symmetric"):`

# Calls

- [SetupEnsightOutputProcess](../../../functions/tests/test_ensight_output_process/SetupEnsightOutputProcess.md)
- [SetSolutionWithMatrixVariables](../../../functions/tests/test_ensight_output_process/SetSolutionWithMatrixVariables.md)
- [CheckFileNonEmpty](../../../functions/tests/test_ensight_output_process/CheckFileNonEmpty.md)

# Called by

- [test_ascii_ensight6_symmetric_tensor_2D](../../../functions/tests/test_ensight_output_process/TestEnsightOutputProcess/test_ascii_ensight6_symmetric_tensor_2D.md)
- [test_binary_ensight6_symmetric_tensor_2D](../../../functions/tests/test_ensight_output_process/TestEnsightOutputProcess/test_binary_ensight6_symmetric_tensor_2D.md)
- [test_ascii_ensightgold_symmetric_tensor_2D](../../../functions/tests/test_ensight_output_process/TestEnsightOutputProcess/test_ascii_ensightgold_symmetric_tensor_2D.md)
- [test_binary_ensightgold_symmetric_tensor_2D](../../../functions/tests/test_ensight_output_process/TestEnsightOutputProcess/test_binary_ensightgold_symmetric_tensor_2D.md)
- [test_ascii_ensight6_symmetric_tensor_3D](../../../functions/tests/test_ensight_output_process/TestEnsightOutputProcess/test_ascii_ensight6_symmetric_tensor_3D.md)
- [test_binary_ensight6_symmetric_tensor_3D](../../../functions/tests/test_ensight_output_process/TestEnsightOutputProcess/test_binary_ensight6_symmetric_tensor_3D.md)
- [test_ascii_ensightgold_symmetric_tensor_3D](../../../functions/tests/test_ensight_output_process/TestEnsightOutputProcess/test_ascii_ensightgold_symmetric_tensor_3D.md)
- [test_binary_ensightgold_symmetric_tensor_3D](../../../functions/tests/test_ensight_output_process/TestEnsightOutputProcess/test_binary_ensightgold_symmetric_tensor_3D.md)
- [test_ascii_ensight6_asymmetric_tensor_2D](../../../functions/tests/test_ensight_output_process/TestEnsightOutputProcess/test_ascii_ensight6_asymmetric_tensor_2D.md)
- [test_binary_ensight6_asymmetric_tensor_2D](../../../functions/tests/test_ensight_output_process/TestEnsightOutputProcess/test_binary_ensight6_asymmetric_tensor_2D.md)
- [test_ascii_ensightgold_asymmetric_tensor_2D](../../../functions/tests/test_ensight_output_process/TestEnsightOutputProcess/test_ascii_ensightgold_asymmetric_tensor_2D.md)
- [test_binary_ensightgold_asymmetric_tensor_2D](../../../functions/tests/test_ensight_output_process/TestEnsightOutputProcess/test_binary_ensightgold_asymmetric_tensor_2D.md)
- [test_ascii_ensight6_asymmetric_tensor_3D](../../../functions/tests/test_ensight_output_process/TestEnsightOutputProcess/test_ascii_ensight6_asymmetric_tensor_3D.md)
- [test_binary_ensight6_asymmetric_tensor_3D](../../../functions/tests/test_ensight_output_process/TestEnsightOutputProcess/test_binary_ensight6_asymmetric_tensor_3D.md)
- [test_ascii_ensightgold_asymmetric_tensor_3D](../../../functions/tests/test_ensight_output_process/TestEnsightOutputProcess/test_ascii_ensightgold_asymmetric_tensor_3D.md)
- [test_binary_ensightgold_asymmetric_tensor_3D](../../../functions/tests/test_ensight_output_process/TestEnsightOutputProcess/test_binary_ensightgold_asymmetric_tensor_3D.md)