---
type: Python Method
title: __define_output_plane
resource: python_scripts/gid_output_process.py#L463-L489
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/GetDouble
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/cutting_utility/CuttingUtility/GenerateCut
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/python_scripts/gid_output_process/GiDOutputProcess/_initialize_cut_output
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def __define_output_plane(self,cut_data):`

# Calls

- [GetDouble](../../../../functions/sources/kratos_parameters/Parameters/GetDouble.md)
- [GenerateCut](../../../../functions/utilities/cutting_utility/CuttingUtility/GenerateCut.md)

# Called by

- [_initialize_cut_output](../../../../functions/python_scripts/gid_output_process/GiDOutputProcess/_initialize_cut_output.md)