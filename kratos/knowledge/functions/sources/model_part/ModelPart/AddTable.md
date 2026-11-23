---
type: C++ Method
title: AddTable
resource: sources/model_part.cpp#L571-L577
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/model_part/IsSubModelPart
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/python_scripts/read_csv_table_utility/ReadCsvTableUtility/Read
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadSubModelPartTablesBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/sources/test_model_part/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_model_part/TestModelPart/test_model_part_tables
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ModelPart::AddTable(ModelPart::IndexType TableId, ModelPart::TableType::Pointer pNewTable)`

# Calls

- [IsSubModelPart](../../../../functions/includes/model_part/IsSubModelPart.md)

# Called by

- [Read](../../../../functions/python_scripts/read_csv_table_utility/ReadCsvTableUtility/Read.md)
- [ReadSubModelPartTablesBlock](../../../../functions/sources/model_part_io/ModelPartIO/ReadSubModelPartTablesBlock.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/sources/test_model_part/KRATOS_TEST_CASE_IN_SUITE.md)
- [test_model_part_tables](../../../../functions/tests/test_model_part/TestModelPart/test_model_part_tables.md)