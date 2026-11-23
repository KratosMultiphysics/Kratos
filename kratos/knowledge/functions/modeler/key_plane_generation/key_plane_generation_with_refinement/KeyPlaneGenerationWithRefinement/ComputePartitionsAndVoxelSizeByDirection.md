---
type: C++ Method
title: ComputePartitionsAndVoxelSizeByDirection
resource: modeler/key_plane_generation/key_plane_generation_with_refinement.cpp#L120-L193
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/modeler/key_plane_generation/key_plane_generation_with_refinement/KeyPlaneGenerationWithRefinement/ReturnLocalTheoreticalVoxelSize
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/modeler/key_plane_generation/key_plane_generation_with_refinement/KeyPlaneGenerationWithRefinement/Generate
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`std::pair<array_1d<std::vector<double>,3>,array_1d<std::vector<double>,3>> KeyPlaneGenerationWithRefinement::ComputePartitionsAndVoxelSizeByDirection(std::vector<BoundingBoxType>& rRefinementAreas, std::vector<array_1d<double,3>>& rVoxelSizes, BoundingBoxType& rGlobalBoundingBox, const array_1d<double,3>& rGlobalVoxelSize)`

# Calls

- [ReturnLocalTheoreticalVoxelSize](../../../../../functions/modeler/key_plane_generation/key_plane_generation_with_refinement/KeyPlaneGenerationWithRefinement/ReturnLocalTheoreticalVoxelSize.md)

# Called by

- [Generate](../../../../../functions/modeler/key_plane_generation/key_plane_generation_with_refinement/KeyPlaneGenerationWithRefinement/Generate.md)