---
type: Python Method
title: _initialize_cut_output
resource: python_scripts/gid_output_process.py#L381-L396
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/cutting_utility/CuttingUtility/FindSmallestEdge
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/cutting_utility/CuttingUtility/AddVariablesToCutModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/ValidateAndAssignDefaults
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/gid_output_process/GiDOutputProcess/__define_output_plane
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/python_scripts/gid_output_process/GiDOutputProcess/ExecuteInitialize
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def _initialize_cut_output(self,plane_output_configuration):`

# Calls

- [FindSmallestEdge](../../../../functions/utilities/cutting_utility/CuttingUtility/FindSmallestEdge.md)
- [AddVariablesToCutModelPart](../../../../functions/utilities/cutting_utility/CuttingUtility/AddVariablesToCutModelPart.md)
- [ValidateAndAssignDefaults](../../../../functions/sources/kratos_parameters/Parameters/ValidateAndAssignDefaults.md)
- [__define_output_plane](../../../../functions/python_scripts/gid_output_process/GiDOutputProcess/__define_output_plane.md)

# Called by

- [ExecuteInitialize](../../../../functions/python_scripts/gid_output_process/GiDOutputProcess/ExecuteInitialize.md)