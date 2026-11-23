---
type: Python Method
title: write_results
resource: python_scripts/gid_output.py#L190-L214
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python_scripts/gid_output/GiDOutput/get_out_model_part
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/cutting_utility/CuttingUtility/UpdateCutData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/gid_output/GiDOutput/_write_mesh
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/gid_output/GiDOutput/_initialize_results
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_interface/kratos_globals/KratosGlobalsImpl/GetVariable
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/gid_output/GiDOutput/_write_nodal_results
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/gid_output/GiDOutput/_write_gp_results
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/gid_output/GiDOutput/_finalize_results
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/gid_output/GiDOutput/write_step_to_list
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def write_results(self, label, model_part, nodal_variables, gp_variables): # label = str(label) #it should be a C double`

# Calls

- [get_out_model_part](../../../../functions/python_scripts/gid_output/GiDOutput/get_out_model_part.md)
- [UpdateCutData](../../../../functions/utilities/cutting_utility/CuttingUtility/UpdateCutData.md)
- [_write_mesh](../../../../functions/python_scripts/gid_output/GiDOutput/_write_mesh.md)
- [_initialize_results](../../../../functions/python_scripts/gid_output/GiDOutput/_initialize_results.md)
- [GetVariable](../../../../functions/python_interface/kratos_globals/KratosGlobalsImpl/GetVariable.md)
- [_write_nodal_results](../../../../functions/python_scripts/gid_output/GiDOutput/_write_nodal_results.md)
- [_write_gp_results](../../../../functions/python_scripts/gid_output/GiDOutput/_write_gp_results.md)
- [_finalize_results](../../../../functions/python_scripts/gid_output/GiDOutput/_finalize_results.md)
- [write_step_to_list](../../../../functions/python_scripts/gid_output/GiDOutput/write_step_to_list.md)