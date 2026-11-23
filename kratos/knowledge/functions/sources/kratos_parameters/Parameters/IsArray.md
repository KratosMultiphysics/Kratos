---
type: C++ Method
title: IsArray
resource: sources/kratos_parameters.cpp#L576-L579
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/includes/json_io/KratosJsonIO/FillModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/json_io/KratosJsonIO/EnforceConsecutiveOrdering
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/cad_json_input/CadJsonInput/ReadBrepSurfaces
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/cad_json_input/CadJsonInput/ReadBrepCurveOnSurfaces
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/cad_json_input/CadJsonInput/ReadPointsOnGeometries
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/cad_json_input/CadJsonInput/ReadControlPointVector
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/serial_model_part_combinator_modeler/SerialModelPartCombinatorModeler/SetupModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/apply_constant_vectorvalue_process/ApplyConstantVectorValueProcess/ApplyConstantVectorValueProcess
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/assign_vector_by_direction_process/AssignVectorByDirectionProcess/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/assign_vector_by_direction_to_entity_process/AssignVectorByDirectionToEntityProcess/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/json_output_process/JsonOutputProcess/__generate_variable_list_from_input
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`bool Parameters::IsArray() const`

# Called by

- [FillModelPart](../../../../functions/includes/json_io/KratosJsonIO/FillModelPart.md)
- [EnforceConsecutiveOrdering](../../../../functions/includes/json_io/KratosJsonIO/EnforceConsecutiveOrdering.md)
- [ReadBrepSurfaces](../../../../functions/input_output/cad_json_input/CadJsonInput/ReadBrepSurfaces.md)
- [ReadBrepCurveOnSurfaces](../../../../functions/input_output/cad_json_input/CadJsonInput/ReadBrepCurveOnSurfaces.md)
- [ReadPointsOnGeometries](../../../../functions/input_output/cad_json_input/CadJsonInput/ReadPointsOnGeometries.md)
- [ReadControlPointVector](../../../../functions/input_output/cad_json_input/CadJsonInput/ReadControlPointVector.md)
- [SetupModelPart](../../../../functions/modeler/serial_model_part_combinator_modeler/SerialModelPartCombinatorModeler/SetupModelPart.md)
- [ApplyConstantVectorValueProcess](../../../../functions/processes/apply_constant_vectorvalue_process/ApplyConstantVectorValueProcess/ApplyConstantVectorValueProcess.md)
- [__init__](../../../../functions/python_scripts/assign_vector_by_direction_process/AssignVectorByDirectionProcess/__init__.md)
- [__init__](../../../../functions/python_scripts/assign_vector_by_direction_to_entity_process/AssignVectorByDirectionToEntityProcess/__init__.md)
- [__generate_variable_list_from_input](../../../../functions/python_scripts/json_output_process/JsonOutputProcess/__generate_variable_list_from_input.md)