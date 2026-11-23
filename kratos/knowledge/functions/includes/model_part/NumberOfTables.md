---
type: C++ Function
title: NumberOfTables
resource: includes/model_part.h#L570-L573
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/sources/model_part/ModelPart/PrintData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_import_model_part_modeler/TestImportMDPAModeler/testImportMDPAModeler
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_model_part/TestModelPart/test_model_part_tables
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_model_part_io/TestModelPartIO/test_model_part_io_read_using_submodelpart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_model_part_io/TestModelPartIO/execute_test_model_part_io_read_model_part
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_model_part_io/TestModelPartIO/test_model_part_io_read_model_part_mesh_only
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_restart/TestRestart/_check_modelpart
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`SizeType NumberOfTables() const`

# Called by

- [PrintData](../../../functions/sources/model_part/ModelPart/PrintData.md)
- [testImportMDPAModeler](../../../functions/tests/test_import_model_part_modeler/TestImportMDPAModeler/testImportMDPAModeler.md)
- [test_model_part_tables](../../../functions/tests/test_model_part/TestModelPart/test_model_part_tables.md)
- [test_model_part_io_read_using_submodelpart](../../../functions/tests/test_model_part_io/TestModelPartIO/test_model_part_io_read_using_submodelpart.md)
- [execute_test_model_part_io_read_model_part](../../../functions/tests/test_model_part_io/TestModelPartIO/execute_test_model_part_io_read_model_part.md)
- [test_model_part_io_read_model_part_mesh_only](../../../functions/tests/test_model_part_io/TestModelPartIO/test_model_part_io_read_model_part_mesh_only.md)
- [_check_modelpart](../../../functions/tests/test_restart/TestRestart/_check_modelpart.md)