---
type: C++ Method
title: PopulateLocalMesh
resource: modeler/combine_model_part_modeler.cpp#L281-L353
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/modeler/combine_model_part_modeler/CombineModelPartModeler/PopulateCommunicators
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/combine_model_part_modeler/CombineModelPartModeler/DuplicateCommunicatorData
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void CombineModelPartModeler::PopulateLocalMesh( Communicator& rReferenceComm, Communicator& rDestinationComm, ModelPart& rDestinationModelPart ) const`

# Called by

- [PopulateCommunicators](../../../../functions/modeler/combine_model_part_modeler/CombineModelPartModeler/PopulateCommunicators.md)
- [DuplicateCommunicatorData](../../../../functions/modeler/combine_model_part_modeler/CombineModelPartModeler/DuplicateCommunicatorData.md)