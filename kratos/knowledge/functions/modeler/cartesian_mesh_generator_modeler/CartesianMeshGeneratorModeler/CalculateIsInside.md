---
type: C++ Method
title: CalculateIsInside
resource: modeler/cartesian_mesh_generator_modeler.h#L306-L335
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/amatrix_interface/Matrix/fill
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/cartesian_mesh_generator_modeler/CartesianMeshGeneratorModeler/NodeIndex
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/modeler/cartesian_mesh_generator_modeler/CartesianMeshGeneratorModeler/GenerateMesh
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void CalculateIsInside(ModelPart& rThisModelPart)`

# Calls

- [fill](../../../../functions/includes/amatrix_interface/Matrix/fill.md)
- [NodeIndex](../../../../functions/modeler/cartesian_mesh_generator_modeler/CartesianMeshGeneratorModeler/NodeIndex.md)

# Called by

- [GenerateMesh](../../../../functions/modeler/cartesian_mesh_generator_modeler/CartesianMeshGeneratorModeler/GenerateMesh.md)