---
type: C++ Method
title: GenerateMesh
resource: modeler/cartesian_mesh_generator_modeler.h#L96-L222
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/point/Point/Dimension
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/timer/Timer/Start
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/cartesian_mesh_generator_modeler/CartesianMeshGeneratorModeler/CalculateDivisionNumbers
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/cartesian_mesh_generator_modeler/CartesianMeshGeneratorModeler/CalculateBoundaryIntersections
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/cartesian_mesh_generator_modeler/CartesianMeshGeneratorModeler/CalculateIsInside
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/cartesian_mesh_generator_modeler/CartesianMeshGeneratorModeler/CalculateNormals
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/cartesian_mesh_generator_modeler/CartesianMeshGeneratorModeler/NodeIndex
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/timer/Timer/Stop
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/cartesian_mesh_generator_modeler/CartesianMeshGeneratorModeler/FindNearestNodeIndex
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void GenerateMesh(ModelPart& rThisModelPart, Element const& rReferenceElement)`

# Calls

- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [Dimension](../../../../functions/geometries/point/Point/Dimension.md)
- [Start](../../../../functions/utilities/timer/Timer/Start.md)
- [CalculateDivisionNumbers](../../../../functions/modeler/cartesian_mesh_generator_modeler/CartesianMeshGeneratorModeler/CalculateDivisionNumbers.md)
- [CalculateBoundaryIntersections](../../../../functions/modeler/cartesian_mesh_generator_modeler/CartesianMeshGeneratorModeler/CalculateBoundaryIntersections.md)
- [CalculateIsInside](../../../../functions/modeler/cartesian_mesh_generator_modeler/CartesianMeshGeneratorModeler/CalculateIsInside.md)
- [CalculateNormals](../../../../functions/modeler/cartesian_mesh_generator_modeler/CartesianMeshGeneratorModeler/CalculateNormals.md)
- [NodeIndex](../../../../functions/modeler/cartesian_mesh_generator_modeler/CartesianMeshGeneratorModeler/NodeIndex.md)
- [Stop](../../../../functions/utilities/timer/Timer/Stop.md)
- [FindNearestNodeIndex](../../../../functions/modeler/cartesian_mesh_generator_modeler/CartesianMeshGeneratorModeler/FindNearestNodeIndex.md)