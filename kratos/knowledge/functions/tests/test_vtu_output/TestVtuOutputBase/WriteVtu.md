---
type: Python Method
title: WriteVtu
resource: tests/test_vtu_output.py#L42-L69
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/vtu_output/VtuOutput/VtuOutput
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/vtu_output/VtuOutput/AddVariable
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/test_vtu_output/TestVtuOutputBase/test_WriteMeshAscii
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_vtu_output/TestVtuOutputBase/test_WriteMeshBinary
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_vtu_output/TestVtuOutputBase/test_WriteMeshRaw
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_vtu_output/TestVtuOutputBase/test_WriteMeshCompressedRaw
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_vtu_output/TestVtuOutputBase/test_WriteMeshAsciiWithError
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_vtu_output/TestVtuOutputBase/test_WriteMeshBinaryWithError
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_vtu_output/TestVtuOutputBase/test_WriteMeshRawWithError
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_vtu_output/TestVtuOutputBase/test_WriteMeshCompressedRawWithError
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def WriteVtu(self, output_format: Kratos.VtuOutput.WriterFormat, error_check = False):`

# Calls

- [VtuOutput](../../../../functions/input_output/vtu_output/VtuOutput/VtuOutput.md)
- [AddVariable](../../../../functions/input_output/vtu_output/VtuOutput/AddVariable.md)

# Called by

- [test_WriteMeshAscii](../../../../functions/tests/test_vtu_output/TestVtuOutputBase/test_WriteMeshAscii.md)
- [test_WriteMeshBinary](../../../../functions/tests/test_vtu_output/TestVtuOutputBase/test_WriteMeshBinary.md)
- [test_WriteMeshRaw](../../../../functions/tests/test_vtu_output/TestVtuOutputBase/test_WriteMeshRaw.md)
- [test_WriteMeshCompressedRaw](../../../../functions/tests/test_vtu_output/TestVtuOutputBase/test_WriteMeshCompressedRaw.md)
- [test_WriteMeshAsciiWithError](../../../../functions/tests/test_vtu_output/TestVtuOutputBase/test_WriteMeshAsciiWithError.md)
- [test_WriteMeshBinaryWithError](../../../../functions/tests/test_vtu_output/TestVtuOutputBase/test_WriteMeshBinaryWithError.md)
- [test_WriteMeshRawWithError](../../../../functions/tests/test_vtu_output/TestVtuOutputBase/test_WriteMeshRawWithError.md)
- [test_WriteMeshCompressedRawWithError](../../../../functions/tests/test_vtu_output/TestVtuOutputBase/test_WriteMeshCompressedRawWithError.md)