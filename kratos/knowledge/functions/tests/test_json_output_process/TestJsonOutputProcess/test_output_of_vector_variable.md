---
type: Python Method
title: test_output_of_vector_variable
resource: tests/test_json_output_process.py#L47-L79
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/tests/test_json_output_process/TestJsonOutputProcess/add_three_nodes_to_model_part
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_json_output_process/TestJsonOutputProcess/add_vectors_to_nodes_of_model_part
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_json_output_process/TestJsonOutputProcess/produce_json_output_file
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/KratosUnittest/TestCase/assertVectorAlmostEqual
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def test_output_of_vector_variable(self):`

# Calls

- [add_three_nodes_to_model_part](../../../../functions/tests/test_json_output_process/TestJsonOutputProcess/add_three_nodes_to_model_part.md)
- [add_vectors_to_nodes_of_model_part](../../../../functions/tests/test_json_output_process/TestJsonOutputProcess/add_vectors_to_nodes_of_model_part.md)
- [produce_json_output_file](../../../../functions/tests/test_json_output_process/TestJsonOutputProcess/produce_json_output_file.md)
- [assertVectorAlmostEqual](../../../../functions/python_scripts/KratosUnittest/TestCase/assertVectorAlmostEqual.md)