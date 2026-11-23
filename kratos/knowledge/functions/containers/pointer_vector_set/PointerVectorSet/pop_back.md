---
type: C++ Method
title: pop_back
resource: containers/pointer_vector_set.h#L556-L561
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/modeler/coloring/color_connected_cells_in_touch/ColorConnectedCellsInTouch/ColorConnectedCellsToThisCell
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/SolveIncludes
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/spatial_containers/cell/Cell/Remove
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/spatial_containers/octree_binary/OctreeBinary/GetLeavesInBoundingBoxNormalized
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/spatial_containers/octree_binary/OctreeBinary/GetAllLeavesVector
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/spatial_containers/octree_binary/OctreeBinary/RefineWithUniformSizeNormalized
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/spatial_containers/octree_binary/OctreeBinary/InsertNormalized
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/spatial_containers/octree_binary/OctreeBinary/Insert
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/spatial_containers/octree_binary/OctreeBinary/GetIntersectedLeaves
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/polynomial_utilities/Divide
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/polynomial_utilities/IsolateRoots
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/string_utilities/StringToVector
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/tessellation_utilities/curve_tessellation/CurveTessellation/ComputeTessellation
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/timer/Timer/Stop
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void pop_back()`

# Called by

- [ColorConnectedCellsToThisCell](../../../../functions/modeler/coloring/color_connected_cells_in_touch/ColorConnectedCellsInTouch/ColorConnectedCellsToThisCell.md)
- [SolveIncludes](../../../../functions/sources/kratos_parameters/Parameters/SolveIncludes.md)
- [Remove](../../../../functions/spatial_containers/cell/Cell/Remove.md)
- [GetLeavesInBoundingBoxNormalized](../../../../functions/spatial_containers/octree_binary/OctreeBinary/GetLeavesInBoundingBoxNormalized.md)
- [GetAllLeavesVector](../../../../functions/spatial_containers/octree_binary/OctreeBinary/GetAllLeavesVector.md)
- [RefineWithUniformSizeNormalized](../../../../functions/spatial_containers/octree_binary/OctreeBinary/RefineWithUniformSizeNormalized.md)
- [InsertNormalized](../../../../functions/spatial_containers/octree_binary/OctreeBinary/InsertNormalized.md)
- [Insert](../../../../functions/spatial_containers/octree_binary/OctreeBinary/Insert.md)
- [GetIntersectedLeaves](../../../../functions/spatial_containers/octree_binary/OctreeBinary/GetIntersectedLeaves.md)
- [Divide](../../../../functions/utilities/polynomial_utilities/Divide.md)
- [IsolateRoots](../../../../functions/utilities/polynomial_utilities/IsolateRoots.md)
- [StringToVector](../../../../functions/utilities/string_utilities/StringToVector.md)
- [ComputeTessellation](../../../../functions/utilities/tessellation_utilities/curve_tessellation/CurveTessellation/ComputeTessellation.md)
- [Stop](../../../../functions/utilities/timer/Timer/Stop.md)