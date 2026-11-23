---
type: Python Function
title: IsMPIAvailable
resource: python_scripts/kratos_utilities.py#L94-L99
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python_scripts/kratos_utilities/GetKratosMultiphysicsPath
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/python_interface/__ModuleInitDetail
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def IsMPIAvailable():`

# Calls

- [GetKratosMultiphysicsPath](../../../functions/python_scripts/kratos_utilities/GetKratosMultiphysicsPath.md)

# Called by

- [__ModuleInitDetail](../../../functions/python_interface/__ModuleInitDetail.md)