---
type: C++ Method
title: CalculateNodalAreaProcess
resource: processes/calculate_nodal_area_process.cpp#L28-L46
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/GetString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/RecursivelyValidateAndAssignDefaults
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetInt
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_nodal_area_process/CalculateNodalAreaProcess<THistorical>/CheckDomainSize
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`CalculateNodalAreaProcess<THistorical>::CalculateNodalAreaProcess( Model& rModel, Parameters ThisParameters ) : mrModelPart(rModel.GetModelPart(ThisParameters["model_part_name"].GetString()))`

# Calls

- [GetString](../../../../functions/sources/kratos_parameters/Parameters/GetString.md)
- [RecursivelyValidateAndAssignDefaults](../../../../functions/sources/kratos_parameters/Parameters/RecursivelyValidateAndAssignDefaults.md)
- [GetInt](../../../../functions/sources/kratos_parameters/Parameters/GetInt.md)
- [CheckDomainSize](../../../../functions/processes/calculate_nodal_area_process/CalculateNodalAreaProcess<THistorical>/CheckDomainSize.md)