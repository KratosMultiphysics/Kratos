---
type: C++ Method
title: GenerateKeyplanes
resource: modeler/key_plane_generation/key_plane_generation_with_refinement.cpp#L213-L234
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/modeler/key_plane_generation/key_plane_generation/VoxelMesherKeyPlaneGeneration/AddKeyPlane
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/modeler/key_plane_generation/key_plane_generation_with_refinement/KeyPlaneGenerationWithRefinement/Generate
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void KeyPlaneGenerationWithRefinement::GenerateKeyplanes(array_1d<std::vector<double>,3>& rPartitionLimits, array_1d<std::vector<double>,3>& rTheoreticalVoxelSize)`

# Calls

- [AddKeyPlane](../../../../../functions/modeler/key_plane_generation/key_plane_generation/VoxelMesherKeyPlaneGeneration/AddKeyPlane.md)

# Called by

- [Generate](../../../../../functions/modeler/key_plane_generation/key_plane_generation_with_refinement/KeyPlaneGenerationWithRefinement/Generate.md)