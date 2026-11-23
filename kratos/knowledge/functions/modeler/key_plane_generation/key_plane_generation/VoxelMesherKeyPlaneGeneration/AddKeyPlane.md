---
type: C++ Method
title: AddKeyPlane
resource: modeler/key_plane_generation/key_plane_generation.cpp#L50-L53
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/modeler/key_plane_generation/key_plane_generation_by_bounding_box/KeyPlaneGenerationByBoundingBox/Generate
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/key_plane_generation/key_plane_generation_with_refinement/KeyPlaneGenerationWithRefinement/GenerateKeyplanes
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void VoxelMesherKeyPlaneGeneration::AddKeyPlane(std::size_t Direction, double Coordinate)`

# Called by

- [Generate](../../../../../functions/modeler/key_plane_generation/key_plane_generation_by_bounding_box/KeyPlaneGenerationByBoundingBox/Generate.md)
- [GenerateKeyplanes](../../../../../functions/modeler/key_plane_generation/key_plane_generation_with_refinement/KeyPlaneGenerationWithRefinement/GenerateKeyplanes.md)