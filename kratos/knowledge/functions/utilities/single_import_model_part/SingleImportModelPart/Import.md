---
type: C++ Method
title: Import
resource: utilities/single_import_model_part.cpp#L26-L75
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/AddMissingParameters
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetBool
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/reorder_consecutive_model_part_io/ReorderConsecutiveModelPartIO
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/modeler/serial_model_part_combinator_modeler/SerialModelPartCombinatorModeler/SetupModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/modelers/import_mdpa_modeler/ImportMDPAModeler/SetupGeometryModel
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/python_solver/PythonSolver/_ImportModelPart
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void SingleImportModelPart::Import( ModelPart& rModelPart, Parameters ModelPartImportParameters, const std::string& InputType )`

# Calls

- [AddMissingParameters](../../../../functions/sources/kratos_parameters/Parameters/AddMissingParameters.md)
- [GetString](../../../../functions/sources/kratos_parameters/Parameters/GetString.md)
- [GetBool](../../../../functions/sources/kratos_parameters/Parameters/GetBool.md)
- [ReorderConsecutiveModelPartIO](../../../../functions/includes/reorder_consecutive_model_part_io/ReorderConsecutiveModelPartIO.md)

# Called by

- [SetupModelPart](../../../../functions/modeler/serial_model_part_combinator_modeler/SerialModelPartCombinatorModeler/SetupModelPart.md)
- [SetupGeometryModel](../../../../functions/python_scripts/modelers/import_mdpa_modeler/ImportMDPAModeler/SetupGeometryModel.md)
- [_ImportModelPart](../../../../functions/python_scripts/python_solver/PythonSolver/_ImportModelPart.md)