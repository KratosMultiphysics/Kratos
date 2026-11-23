---
type: C++ Method
title: AddFlag
resource: input_output/vtu_output.cpp#L823-L865
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/python_scripts/vtu_output_process/VtuOutputProcess/ExecuteInitialize
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void VtuOutput::AddFlag( const std::string& rFlagName, const Flags& rFlagVariable, const Globals::DataLocation& DataLocation)`

# Called by

- [ExecuteInitialize](../../../../functions/python_scripts/vtu_output_process/VtuOutputProcess/ExecuteInitialize.md)