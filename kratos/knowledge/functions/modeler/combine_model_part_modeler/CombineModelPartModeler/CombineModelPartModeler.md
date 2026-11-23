---
type: C++ Method
title: CombineModelPartModeler
resource: modeler/combine_model_part_modeler.cpp#L16-L24
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/ValidateAndAssignDefaults
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`CombineModelPartModeler::CombineModelPartModeler( Model& rModel, Parameters ModelerParameters ) : Modeler(rModel, ModelerParameters), mpModel(&rModel), mParameters(ModelerParameters)`

# Calls

- [ValidateAndAssignDefaults](../../../../functions/sources/kratos_parameters/Parameters/ValidateAndAssignDefaults.md)