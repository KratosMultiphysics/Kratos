---
type: Python Method
title: test_CombineModelParts
resource: tests/test_combine_model_part_modeler.py#L53-L164
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/tests/test_combine_model_part_modeler/TestCombineModelPartModeler/_CreateModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/SetString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/project/Project/GetModel
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/communicator/Communicator/GlobalNumberOfNodes
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/communicator/Communicator/GlobalNumberOfElements
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/communicator/Communicator/GlobalNumberOfConditions
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def test_CombineModelParts(self):`

# Calls

- [_CreateModelPart](../../../../functions/tests/test_combine_model_part_modeler/TestCombineModelPartModeler/_CreateModelPart.md)
- [GetString](../../../../functions/sources/kratos_parameters/Parameters/GetString.md)
- [SetString](../../../../functions/sources/kratos_parameters/Parameters/SetString.md)
- [GetModel](../../../../functions/python_scripts/project/Project/GetModel.md)
- [GlobalNumberOfNodes](../../../../functions/sources/communicator/Communicator/GlobalNumberOfNodes.md)
- [GlobalNumberOfElements](../../../../functions/sources/communicator/Communicator/GlobalNumberOfElements.md)
- [GlobalNumberOfConditions](../../../../functions/sources/communicator/Communicator/GlobalNumberOfConditions.md)