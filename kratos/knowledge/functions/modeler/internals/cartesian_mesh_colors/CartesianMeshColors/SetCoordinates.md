---
type: C++ Method
title: SetCoordinates
resource: modeler/internals/cartesian_mesh_colors.h#L69-L102
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/modeler/internals/cartesian_mesh_colors/CartesianMeshColors/SetAllColors
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/modeler/voxel_mesh_generator_modeler/VoxelMeshGeneratorModeler/PreparingTheInternalDataStructure
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/modeler/test_voxel_mesh_modeler/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void SetCoordinates(std::vector<double> const& rNodalXCoordinates, std::vector<double> const& rNodalYCoordinates, std::vector<double> const& rNodalZCoordinates)`

# Calls

- [SetAllColors](../../../../../functions/modeler/internals/cartesian_mesh_colors/CartesianMeshColors/SetAllColors.md)

# Called by

- [PreparingTheInternalDataStructure](../../../../../functions/modeler/voxel_mesh_generator_modeler/VoxelMeshGeneratorModeler/PreparingTheInternalDataStructure.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../../functions/tests/cpp_tests/modeler/test_voxel_mesh_modeler/KRATOS_TEST_CASE_IN_SUITE.md)