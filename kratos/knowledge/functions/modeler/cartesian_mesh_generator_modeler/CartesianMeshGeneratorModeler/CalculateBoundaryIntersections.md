---
type: C++ Method
title: CalculateBoundaryIntersections
resource: modeler/cartesian_mesh_generator_modeler.h#L337-L388
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/modeler/cartesian_mesh_generator_modeler/CartesianMeshGeneratorModeler/GenerateMesh
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void CalculateBoundaryIntersections(ModelPart& rThisModelPart)`

# Calls

- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)

# Called by

- [GenerateMesh](../../../../functions/modeler/cartesian_mesh_generator_modeler/CartesianMeshGeneratorModeler/GenerateMesh.md)