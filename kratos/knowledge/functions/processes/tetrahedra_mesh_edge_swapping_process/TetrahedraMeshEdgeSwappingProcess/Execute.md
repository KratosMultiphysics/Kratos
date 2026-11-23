---
type: C++ Method
title: Execute
resource: processes/tetrahedra_mesh_edge_swapping_process.cpp#L38-L72
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/tetrahedra_edge_shell/TetrahedraEdgeShell/AddTetrahedron
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/tetrahedra_edge_shell/TetrahedraEdgeShell/GetNumberOfTetrahedra
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void TetrahedraMeshEdgeSwappingProcess::Execute()`

# Calls

- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [AddTetrahedron](../../../../functions/modeler/tetrahedra_edge_shell/TetrahedraEdgeShell/AddTetrahedron.md)
- [GetNumberOfTetrahedra](../../../../functions/modeler/tetrahedra_edge_shell/TetrahedraEdgeShell/GetNumberOfTetrahedra.md)