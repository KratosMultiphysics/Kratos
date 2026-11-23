---
type: C++ Method
title: AddVariablesToCutModelPart
resource: utilities/cutting_utility.h#L312-L319
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/python_scripts/gid_output/GiDOutput/define_cuts
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/gid_output_process/GiDOutputProcess/_initialize_cut_output
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void AddVariablesToCutModelPart( const ModelPart& rModelPart, ModelPart& rNewModelPart) const`

# Called by

- [define_cuts](../../../../functions/python_scripts/gid_output/GiDOutput/define_cuts.md)
- [_initialize_cut_output](../../../../functions/python_scripts/gid_output_process/GiDOutputProcess/_initialize_cut_output.md)