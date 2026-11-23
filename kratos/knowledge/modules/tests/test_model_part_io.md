---
type: Python Module
title: test_model_part_io
resource: tests/test_model_part_io.py#L1-L411
generated:
  by: okf-rs/0.7.0
relationships:
  imports:
  - target: external/os
    resolved_by: tree-sitter
    confidence: exact
  - target: external/pathlib
    resolved_by: tree-sitter
    confidence: exact
  - target: external/kratosmultiphysics
    resolved_by: tree-sitter
    confidence: exact
  - target: external/kratosmultiphysics-kratosunittest-as-kratosunittest
    resolved_by: tree-sitter
    confidence: exact
  - target: external/kratosmultiphysics-kratos-utilities-as-kratosutils
    resolved_by: tree-sitter
    confidence: exact
  - target: external/kratosmultiphysics-testing-utilities
    resolved_by: tree-sitter
    confidence: exact
  - target: external/kratosmultiphysics-structuralmechanicsapplication
    resolved_by: tree-sitter
    confidence: exact
  - target: external/filecmp
    resolved_by: tree-sitter
    confidence: exact
  member_of:
  - target: packages/KratosMultiphysics
    resolved_by: tree-sitter
    confidence: exact
---

# Contains

- [GetFilePath](../../functions/tests/test_model_part_io/GetFilePath.md)
- [TestModelPartIO](../../classes/tests/test_model_part_io/TestModelPartIO.md)
- [tearDown](../../functions/tests/test_model_part_io/TestModelPartIO/tearDown.md)
- [test_model_part_io_read_model_part](../../functions/tests/test_model_part_io/TestModelPartIO/test_model_part_io_read_model_part.md)
- [test_model_part_io_read_using_submodelpart](../../functions/tests/test_model_part_io/TestModelPartIO/test_model_part_io_read_using_submodelpart.md)
- [execute_test_model_part_io_read_model_part](../../functions/tests/test_model_part_io/TestModelPartIO/execute_test_model_part_io_read_model_part.md)
- [test_model_part_io_read_model_part_mesh_only](../../functions/tests/test_model_part_io/TestModelPartIO/test_model_part_io_read_model_part_mesh_only.md)
- [test_model_part_io_write_model_part](../../functions/tests/test_model_part_io/TestModelPartIO/test_model_part_io_write_model_part.md)
- [test_model_part_io_write_model_part_mesh_only](../../functions/tests/test_model_part_io/TestModelPartIO/test_model_part_io_write_model_part_mesh_only.md)
- [test_error_on_wrong_input](../../functions/tests/test_model_part_io/TestModelPartIO/test_error_on_wrong_input.md)
- [TestModelPartIOMPI](../../classes/tests/test_model_part_io/TestModelPartIOMPI.md)
- [test_model_part_io_read_entity_data](../../functions/tests/test_model_part_io/TestModelPartIOMPI/test_model_part_io_read_entity_data.md)
- [GetScalar](../../functions/tests/test_model_part_io/GetScalar.md)
- [GetArray3](../../functions/tests/test_model_part_io/GetArray3.md)
- [GetVector](../../functions/tests/test_model_part_io/GetVector.md)
- [GetMatrix](../../functions/tests/test_model_part_io/GetMatrix.md)
- [CheckEntities](../../functions/tests/test_model_part_io/CheckEntities.md)

# Imports

- `os`
- `pathlib`
- `KratosMultiphysics`
- `KratosMultiphysics.KratosUnittest as KratosUnittest`
- `KratosMultiphysics.kratos_utilities as KratosUtils`
- `KratosMultiphysics.testing.utilities`
- `KratosMultiphysics.StructuralMechanicsApplication`
- `filecmp`

# Member of

- [KratosMultiphysics](../../packages/KratosMultiphysics.md)