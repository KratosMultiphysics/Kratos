---
type: C++ Method
title: NodeIndex
resource: modeler/cartesian_mesh_generator_modeler.h#L459-L462
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/modeler/cartesian_mesh_generator_modeler/CartesianMeshGeneratorModeler/GenerateMesh
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/cartesian_mesh_generator_modeler/CartesianMeshGeneratorModeler/FindNearestNodeIndex
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/cartesian_mesh_generator_modeler/CartesianMeshGeneratorModeler/CalculateIsInside
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`unsigned int NodeIndex(unsigned int i, unsigned int j, unsigned int k)`

# Called by

- [GenerateMesh](../../../../functions/modeler/cartesian_mesh_generator_modeler/CartesianMeshGeneratorModeler/GenerateMesh.md)
- [FindNearestNodeIndex](../../../../functions/modeler/cartesian_mesh_generator_modeler/CartesianMeshGeneratorModeler/FindNearestNodeIndex.md)
- [CalculateIsInside](../../../../functions/modeler/cartesian_mesh_generator_modeler/CartesianMeshGeneratorModeler/CalculateIsInside.md)