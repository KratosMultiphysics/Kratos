---
type: C++ Method
title: VtkOutput
resource: input_output/vtk_output.cpp#L66-L128
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/ValidateAndAssignDefaults
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetInt
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/communicator/Communicator/GlobalNumberOfElements
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/communicator/Communicator/GlobalNumberOfConditions
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/python_scripts/vtk_output_process/VtkOutputProcess/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/exact_mortar_segmentation_utility/ExactMortarIntegrationUtility<TDim, TNumNodes, TBelong, TNumNodesMaster>/TestIODebug
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`VtkOutput::VtkOutput( ModelPart& rModelPart, Parameters ThisParameters ) : mrModelPart(rModelPart), mOutputSettings(ThisParameters)`

# Calls

- [ValidateAndAssignDefaults](../../../../functions/sources/kratos_parameters/Parameters/ValidateAndAssignDefaults.md)
- [GetInt](../../../../functions/sources/kratos_parameters/Parameters/GetInt.md)
- [GetString](../../../../functions/sources/kratos_parameters/Parameters/GetString.md)
- [GlobalNumberOfElements](../../../../functions/sources/communicator/Communicator/GlobalNumberOfElements.md)
- [GlobalNumberOfConditions](../../../../functions/sources/communicator/Communicator/GlobalNumberOfConditions.md)

# Called by

- [__init__](../../../../functions/python_scripts/vtk_output_process/VtkOutputProcess/__init__.md)
- [TestIODebug](../../../../functions/utilities/exact_mortar_segmentation_utility/ExactMortarIntegrationUtility<TDim, TNumNodes, TBelong, TNumNodesMaster>/TestIODebug.md)