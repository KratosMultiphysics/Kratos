---
type: Python Function
title: ReadModelPartFromOBJ
resource: tests/test_obj_io.py#L125-L141
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/SetBool
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/test_obj_io/TestObjIO/test_ReadObjIO
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_obj_io/TestObjIO/test_ReadObjIOTriangles
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_obj_io/TestObjIO/test_ReadObjIOTrianglesDegenerated
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_obj_io/TestObjIO/test_WriteObjIO
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def ReadModelPartFromOBJ(model_part, obj_file, decompose_quad=False):`

# Calls

- [SetBool](../../../functions/sources/kratos_parameters/Parameters/SetBool.md)

# Called by

- [test_ReadObjIO](../../../functions/tests/test_obj_io/TestObjIO/test_ReadObjIO.md)
- [test_ReadObjIOTriangles](../../../functions/tests/test_obj_io/TestObjIO/test_ReadObjIOTriangles.md)
- [test_ReadObjIOTrianglesDegenerated](../../../functions/tests/test_obj_io/TestObjIO/test_ReadObjIOTrianglesDegenerated.md)
- [test_WriteObjIO](../../../functions/tests/test_obj_io/TestObjIO/test_WriteObjIO.md)