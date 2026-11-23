---
type: C++ Method
title: FindSmallestEdge
resource: utilities/cutting_utility.h#L100-L129
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/python_scripts/gid_output/GiDOutput/define_cuts
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/gid_output_process/GiDOutputProcess/_initialize_cut_output
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void FindSmallestEdge(ModelPart& mr_model_part)`

# Calls

- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)

# Called by

- [define_cuts](../../../../functions/python_scripts/gid_output/GiDOutput/define_cuts.md)
- [_initialize_cut_output](../../../../functions/python_scripts/gid_output_process/GiDOutputProcess/_initialize_cut_output.md)