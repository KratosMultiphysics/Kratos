---
type: C++ Method
title: CalculateMortarOperators
resource: includes/mortar_classes.h#L1107-L1127
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/processes/simple_mortar_mapper_process/SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>/AssemblyMortarOperators
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void CalculateMortarOperators( KinematicVariables& rKinematicVariables, const double rIntegrationWeight )`

# Called by

- [AssemblyMortarOperators](../../../../functions/processes/simple_mortar_mapper_process/SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>/AssemblyMortarOperators.md)