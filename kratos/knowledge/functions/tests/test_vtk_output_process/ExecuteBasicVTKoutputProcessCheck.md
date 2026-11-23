---
type: Python Function
title: ExecuteBasicVTKoutputProcessCheck
resource: tests/test_vtk_output_process.py#L407-L476
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/tests/test_vtk_output_process/SetupModelPartHexahedra3D27N
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/SetString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_vtk_output_process/SetupVtkOutputProcess
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/test_vtk_output_process/TestVtkOutputProcess/test_ascii_vtk_output_2D
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_vtk_output_process/TestVtkOutputProcess/test_binary_vtk_output_2D
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_vtk_output_process/TestVtkOutputProcess/test_ascii_vtk_output_3D
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_vtk_output_process/TestVtkOutputProcess/test_binary_vtk_output_3D
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_vtk_output_process/TestVtkOutputProcess/test_ascii_vtk_output_quad_3D
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_vtk_output_process/TestVtkOutputProcess/test_binary_vtk_output_quad_3D
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_vtk_output_process/TestVtkOutputProcess/test_ascii_vtk_output_quad_prism_3D
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_vtk_output_process/TestVtkOutputProcess/test_binary_vtk_output_quad_prism_3D
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_vtk_output_process/TestVtkOutputProcess/test_ascii_vtk_output_quad_hexahedra_3D
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_vtk_output_process/TestVtkOutputProcess/test_binary_vtk_output_quad_hexahedra_3D
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_vtk_output_process/TestVtkOutputProcess/test_ascii_vtk_output_quad_hexahedra_3D_27N
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_vtk_output_process/TestVtkOutputProcess/test_binary_vtk_output_quad_hexahedra_3D_27N
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def ExecuteBasicVTKoutputProcessCheck(file_format = "ascii", setup = "2D"):`

# Calls

- [SetupModelPartHexahedra3D27N](../../../functions/tests/test_vtk_output_process/SetupModelPartHexahedra3D27N.md)
- [SetString](../../../functions/sources/kratos_parameters/Parameters/SetString.md)
- [SetupVtkOutputProcess](../../../functions/tests/test_vtk_output_process/SetupVtkOutputProcess.md)

# Called by

- [test_ascii_vtk_output_2D](../../../functions/tests/test_vtk_output_process/TestVtkOutputProcess/test_ascii_vtk_output_2D.md)
- [test_binary_vtk_output_2D](../../../functions/tests/test_vtk_output_process/TestVtkOutputProcess/test_binary_vtk_output_2D.md)
- [test_ascii_vtk_output_3D](../../../functions/tests/test_vtk_output_process/TestVtkOutputProcess/test_ascii_vtk_output_3D.md)
- [test_binary_vtk_output_3D](../../../functions/tests/test_vtk_output_process/TestVtkOutputProcess/test_binary_vtk_output_3D.md)
- [test_ascii_vtk_output_quad_3D](../../../functions/tests/test_vtk_output_process/TestVtkOutputProcess/test_ascii_vtk_output_quad_3D.md)
- [test_binary_vtk_output_quad_3D](../../../functions/tests/test_vtk_output_process/TestVtkOutputProcess/test_binary_vtk_output_quad_3D.md)
- [test_ascii_vtk_output_quad_prism_3D](../../../functions/tests/test_vtk_output_process/TestVtkOutputProcess/test_ascii_vtk_output_quad_prism_3D.md)
- [test_binary_vtk_output_quad_prism_3D](../../../functions/tests/test_vtk_output_process/TestVtkOutputProcess/test_binary_vtk_output_quad_prism_3D.md)
- [test_ascii_vtk_output_quad_hexahedra_3D](../../../functions/tests/test_vtk_output_process/TestVtkOutputProcess/test_ascii_vtk_output_quad_hexahedra_3D.md)
- [test_binary_vtk_output_quad_hexahedra_3D](../../../functions/tests/test_vtk_output_process/TestVtkOutputProcess/test_binary_vtk_output_quad_hexahedra_3D.md)
- [test_ascii_vtk_output_quad_hexahedra_3D_27N](../../../functions/tests/test_vtk_output_process/TestVtkOutputProcess/test_ascii_vtk_output_quad_hexahedra_3D_27N.md)
- [test_binary_vtk_output_quad_hexahedra_3D_27N](../../../functions/tests/test_vtk_output_process/TestVtkOutputProcess/test_binary_vtk_output_quad_hexahedra_3D_27N.md)