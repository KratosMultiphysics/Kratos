---
type: Python Method
title: _mapper_tests
resource: tests/test_mortar_mapper.py#L80-L133
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/tests/test_mortar_mapper/TestMortarMapperCore/__base_test_mapping
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/SetDouble
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetDouble
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/SetString
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/test_mortar_mapper/TestMortarMapperCore/test_less_basic_mortar_mapping_triangle_pure_implicit
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_mortar_mapper/TestMortarMapperCore/test_less_basic_mortar_mapping_triangle
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_mortar_mapper/TestMortarMapperCore/test_simple_curvature_mortar_mapping_triangle
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_mortar_mapper/TestMortarMapperCore/test_mortar_mapping_triangle
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_mortar_mapper/TestMortarMapperCore/test_mortar_mapping_triangle_discontinous_interface
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_mortar_mapper/TestMortarMapperCore/test_mortar_mapping_quad
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_mortar_mapper/TestMortarMapperCore/test_mortar_mapping_quad_tri
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_mortar_mapper/TestMortarMapperCore/test_mortar_mapping_tri_quad
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def _mapper_tests(self, input_filename, num_nodes, master_num_nodes, pure_implicit = False, inverted = False, discontinuous = False, origin_are_conditions = True, destination_are_conditions = True, consider_tessellation = False, tolerance_factor = 1.0):`

# Calls

- [__base_test_mapping](../../../../functions/tests/test_mortar_mapper/TestMortarMapperCore/__base_test_mapping.md)
- [SetDouble](../../../../functions/sources/kratos_parameters/Parameters/SetDouble.md)
- [GetDouble](../../../../functions/sources/kratos_parameters/Parameters/GetDouble.md)
- [SetString](../../../../functions/sources/kratos_parameters/Parameters/SetString.md)

# Called by

- [test_less_basic_mortar_mapping_triangle_pure_implicit](../../../../functions/tests/test_mortar_mapper/TestMortarMapperCore/test_less_basic_mortar_mapping_triangle_pure_implicit.md)
- [test_less_basic_mortar_mapping_triangle](../../../../functions/tests/test_mortar_mapper/TestMortarMapperCore/test_less_basic_mortar_mapping_triangle.md)
- [test_simple_curvature_mortar_mapping_triangle](../../../../functions/tests/test_mortar_mapper/TestMortarMapperCore/test_simple_curvature_mortar_mapping_triangle.md)
- [test_mortar_mapping_triangle](../../../../functions/tests/test_mortar_mapper/TestMortarMapperCore/test_mortar_mapping_triangle.md)
- [test_mortar_mapping_triangle_discontinous_interface](../../../../functions/tests/test_mortar_mapper/TestMortarMapperCore/test_mortar_mapping_triangle_discontinous_interface.md)
- [test_mortar_mapping_quad](../../../../functions/tests/test_mortar_mapper/TestMortarMapperCore/test_mortar_mapping_quad.md)
- [test_mortar_mapping_quad_tri](../../../../functions/tests/test_mortar_mapper/TestMortarMapperCore/test_mortar_mapping_quad_tri.md)
- [test_mortar_mapping_tri_quad](../../../../functions/tests/test_mortar_mapper/TestMortarMapperCore/test_mortar_mapping_tri_quad.md)