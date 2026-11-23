---
type: C++ Method
title: FindNearestNodeIndex
resource: modeler/cartesian_mesh_generator_modeler.h#L254-L304
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/modeler/cartesian_mesh_generator_modeler/CartesianMeshGeneratorModeler/NodeIndex
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/modeler/cartesian_mesh_generator_modeler/CartesianMeshGeneratorModeler/GenerateMesh
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`unsigned int FindNearestNodeIndex(PointType& rThisPoint, array_1d<double,3>& rNormal)`

# Calls

- [NodeIndex](../../../../functions/modeler/cartesian_mesh_generator_modeler/CartesianMeshGeneratorModeler/NodeIndex.md)

# Called by

- [GenerateMesh](../../../../functions/modeler/cartesian_mesh_generator_modeler/CartesianMeshGeneratorModeler/GenerateMesh.md)