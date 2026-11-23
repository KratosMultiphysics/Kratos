---
type: C++ Method
title: WriteModelPartWithoutNodesToFile
resource: input_output/vtk_output.cpp#L1137-L1179
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python_scripts/project/Project/GetModel
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/fast_transfer_between_model_parts_process/FastTransferBetweenModelPartsProcess/FastTransferBetweenModelPartsProcess
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part/ModelPart/RemoveNodes
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/vtk_output/VtkOutput/WriteModelPartToFile
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/model/Model/DeleteModelPart
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/input_output/vtk_output/VtkOutput/PrintOutput
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void VtkOutput::WriteModelPartWithoutNodesToFile(ModelPart& rModelPart, const std::string& rOutputFilename)`

# Calls

- [GetModel](../../../../functions/python_scripts/project/Project/GetModel.md)
- [FastTransferBetweenModelPartsProcess](../../../../functions/processes/fast_transfer_between_model_parts_process/FastTransferBetweenModelPartsProcess/FastTransferBetweenModelPartsProcess.md)
- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [RemoveNodes](../../../../functions/sources/model_part/ModelPart/RemoveNodes.md)
- [WriteModelPartToFile](../../../../functions/input_output/vtk_output/VtkOutput/WriteModelPartToFile.md)
- [DeleteModelPart](../../../../functions/containers/model/Model/DeleteModelPart.md)

# Called by

- [PrintOutput](../../../../functions/input_output/vtk_output/VtkOutput/PrintOutput.md)