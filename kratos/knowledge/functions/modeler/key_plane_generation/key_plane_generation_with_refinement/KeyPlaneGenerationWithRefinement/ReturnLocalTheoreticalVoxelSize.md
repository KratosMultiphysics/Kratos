---
type: C++ Method
title: ReturnLocalTheoreticalVoxelSize
resource: modeler/key_plane_generation/key_plane_generation_with_refinement.cpp#L195-L211
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/modeler/key_plane_generation/key_plane_generation_with_refinement/KeyPlaneGenerationWithRefinement/ComputePartitionsAndVoxelSizeByDirection
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`double KeyPlaneGenerationWithRefinement::ReturnLocalTheoreticalVoxelSize(std::vector<BoundingBoxType>& rRefinementAreas, std::vector<array_1d<double,3>>& rVoxelSizes, BoundingBoxType& rGlobalBoundingBox, const array_1d<double,3>& rGlobalVoxelSize, int Direction, double Position)`

# Called by

- [ComputePartitionsAndVoxelSizeByDirection](../../../../../functions/modeler/key_plane_generation/key_plane_generation_with_refinement/KeyPlaneGenerationWithRefinement/ComputePartitionsAndVoxelSizeByDirection.md)