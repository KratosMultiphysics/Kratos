---
type: Python Method
title: __base_test_mapping
resource: tests/test_mortar_mapper.py#L18-L78
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/SetBool
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/SetString
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/test_mortar_mapper/TestMortarMapperCore/_mapper_tests
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def __base_test_mapping(self, input_filename, num_nodes, master_num_nodes, pure_implicit, inverted, discontinuous, origin_are_conditions, destination_are_conditions, consider_tessellation):`

# Calls

- [SetBool](../../../../functions/sources/kratos_parameters/Parameters/SetBool.md)
- [SetString](../../../../functions/sources/kratos_parameters/Parameters/SetString.md)

# Called by

- [_mapper_tests](../../../../functions/tests/test_mortar_mapper/TestMortarMapperCore/_mapper_tests.md)