---
type: C++ Method
title: FillAuxiliaryModelPart
resource: processes/skin_detection_process.cpp#L184-L273
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/GetString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetInt
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/data_communicator/Send
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/data_communicator/Recv
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/variable_utils/SetFlag
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/processes/skin_detection_process/SkinDetectionProcess<TDim>/Execute
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void SkinDetectionProcess<TDim>::FillAuxiliaryModelPart( ModelPart& rAuxiliaryModelPart, HashMapVectorIntType& rInverseFaceMap, HashMapVectorIntIdsType& rPropertiesFaceMap )`

# Calls

- [GetString](../../../../functions/sources/kratos_parameters/Parameters/GetString.md)
- [GetInt](../../../../functions/sources/kratos_parameters/Parameters/GetInt.md)
- [Send](../../../../functions/includes/data_communicator/Send.md)
- [Recv](../../../../functions/includes/data_communicator/Recv.md)
- [SetFlag](../../../../functions/utilities/variable_utils/SetFlag.md)

# Called by

- [Execute](../../../../functions/processes/skin_detection_process/SkinDetectionProcess<TDim>/Execute.md)