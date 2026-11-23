---
type: Python Method
title: ConstructItem
resource: python_scripts/model_parameters_factory.py#L16-L81
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/GetString
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/python_scripts/model_parameters_factory/KratosModelParametersFactory/ConstructListOfItems
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/output_process_with_controller/OutputProcessWithController/__init__
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def ConstructItem(self, item): # Check the input settings # Note that we support both 'parameters' and 'Parameters' to keep processes backward compatibility`

# Calls

- [GetString](../../../../functions/sources/kratos_parameters/Parameters/GetString.md)

# Called by

- [ConstructListOfItems](../../../../functions/python_scripts/model_parameters_factory/KratosModelParametersFactory/ConstructListOfItems.md)
- [__init__](../../../../functions/python_scripts/output_process_with_controller/OutputProcessWithController/__init__.md)