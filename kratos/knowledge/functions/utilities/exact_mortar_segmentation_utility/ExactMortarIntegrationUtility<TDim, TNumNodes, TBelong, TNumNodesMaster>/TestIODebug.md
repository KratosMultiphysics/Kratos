---
type: C++ Method
title: TestIODebug
resource: utilities/exact_mortar_segmentation_utility.cpp#L1157-L1318
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python_scripts/project/Project/GetModel
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part/ModelPart/CreateNewProperties
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/flags/IsDefined
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/flags/IsNot
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/gid_io/InitializeMesh
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/gid_io/FinalizeMesh
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/gid_io/CloseResultFile
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/vtk_output/VtkOutput/VtkOutput
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/model/Model/DeleteModelPart
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ExactMortarIntegrationUtility<TDim, TNumNodes, TBelong, TNumNodesMaster>::TestIODebug( ModelPart& rMainModelPart, const std::string IOConsidered )`

# Calls

- [GetModel](../../../../functions/python_scripts/project/Project/GetModel.md)
- [CreateNewProperties](../../../../functions/sources/model_part/ModelPart/CreateNewProperties.md)
- [IsDefined](../../../../functions/containers/flags/IsDefined.md)
- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [IsNot](../../../../functions/containers/flags/IsNot.md)
- [InitializeMesh](../../../../functions/includes/gid_io/InitializeMesh.md)
- [FinalizeMesh](../../../../functions/includes/gid_io/FinalizeMesh.md)
- [CloseResultFile](../../../../functions/includes/gid_io/CloseResultFile.md)
- [VtkOutput](../../../../functions/input_output/vtk_output/VtkOutput/VtkOutput.md)
- [DeleteModelPart](../../../../functions/containers/model/Model/DeleteModelPart.md)