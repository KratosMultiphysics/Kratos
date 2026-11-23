---
type: Python Method
title: __SearchPoint
resource: python_scripts/point_output_process.py#L143-L192
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/GetString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/brute_force_point_locator/BruteForcePointLocator
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/brute_force_point_locator/BruteForcePointLocator/FindNode
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/brute_force_point_locator/BruteForcePointLocator/FindElement
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/brute_force_point_locator/BruteForcePointLocator/FindCondition
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/point_output_process/GetFileHeader
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/python_scripts/point_output_process/PointOutputProcess/ExecuteInitialize
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/point_output_process/PointOutputProcess/ExecuteInitializeSolutionStep
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def __SearchPoint(self): # retrieving the entity type`

# Calls

- [GetString](../../../../functions/sources/kratos_parameters/Parameters/GetString.md)
- [BruteForcePointLocator](../../../../functions/utilities/brute_force_point_locator/BruteForcePointLocator.md)
- [FindNode](../../../../functions/utilities/brute_force_point_locator/BruteForcePointLocator/FindNode.md)
- [FindElement](../../../../functions/utilities/brute_force_point_locator/BruteForcePointLocator/FindElement.md)
- [FindCondition](../../../../functions/utilities/brute_force_point_locator/BruteForcePointLocator/FindCondition.md)
- [GetFileHeader](../../../../functions/python_scripts/point_output_process/GetFileHeader.md)

# Called by

- [ExecuteInitialize](../../../../functions/python_scripts/point_output_process/PointOutputProcess/ExecuteInitialize.md)
- [ExecuteInitializeSolutionStep](../../../../functions/python_scripts/point_output_process/PointOutputProcess/ExecuteInitializeSolutionStep.md)