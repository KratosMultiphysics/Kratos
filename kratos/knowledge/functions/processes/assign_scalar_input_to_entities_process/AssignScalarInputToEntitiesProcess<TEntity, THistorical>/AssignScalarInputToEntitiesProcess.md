---
type: C++ Method
title: AssignScalarInputToEntitiesProcess
resource: processes/assign_scalar_input_to_entities_process.cpp#L36-L41
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/GetString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/ValidateAndAssignDefaults
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/assign_scalar_input_to_entities_process/ConvertAlgorithmString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/string_utilities/ContainsPartialString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/assign_scalar_input_to_entities_process/AssignScalarInputToEntitiesProcess<TEntity, THistorical>/IdentifyDataTXT
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/assign_scalar_input_to_entities_process/AssignScalarInputToEntitiesProcess<TEntity, THistorical>/IdentifyDataJSON
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/assign_scalar_input_to_entities_process/AssignScalarInputToEntitiesProcess<TEntity, THistorical>/ReadDataTXT
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/assign_scalar_input_to_entities_process/AssignScalarInputToEntitiesProcess<TEntity, THistorical>/ReadDataJSON
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/assign_scalar_input_to_entities_process/AssignScalarInputToEntitiesProcess<TEntity, THistorical>/ComputeExtrapolationWeight
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`AssignScalarInputToEntitiesProcess<TEntity, THistorical>::AssignScalarInputToEntitiesProcess( Model& rModel, Parameters rParameters ) : AssignScalarInputToEntitiesProcess(rModel.GetModelPart(rParameters["model_part_name"].GetString()), rParameters)`

# Calls

- [GetString](../../../../functions/sources/kratos_parameters/Parameters/GetString.md)
- [ValidateAndAssignDefaults](../../../../functions/sources/kratos_parameters/Parameters/ValidateAndAssignDefaults.md)
- [ConvertAlgorithmString](../../../../functions/processes/assign_scalar_input_to_entities_process/ConvertAlgorithmString.md)
- [ContainsPartialString](../../../../functions/utilities/string_utilities/ContainsPartialString.md)
- [IdentifyDataTXT](../../../../functions/processes/assign_scalar_input_to_entities_process/AssignScalarInputToEntitiesProcess<TEntity, THistorical>/IdentifyDataTXT.md)
- [IdentifyDataJSON](../../../../functions/processes/assign_scalar_input_to_entities_process/AssignScalarInputToEntitiesProcess<TEntity, THistorical>/IdentifyDataJSON.md)
- [ReadDataTXT](../../../../functions/processes/assign_scalar_input_to_entities_process/AssignScalarInputToEntitiesProcess<TEntity, THistorical>/ReadDataTXT.md)
- [ReadDataJSON](../../../../functions/processes/assign_scalar_input_to_entities_process/AssignScalarInputToEntitiesProcess<TEntity, THistorical>/ReadDataJSON.md)
- [ComputeExtrapolationWeight](../../../../functions/processes/assign_scalar_input_to_entities_process/AssignScalarInputToEntitiesProcess<TEntity, THistorical>/ComputeExtrapolationWeight.md)