---
type: C++ Method
title: IsCompatibleVariable
resource: input_output/vtk_output.cpp#L464-L485
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/input_output/vtk_output/VtkOutput/WriteNodalResultsToFile
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/vtk_output/VtkOutput/WriteElementResultsToFile
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/vtk_output/VtkOutput/WriteConditionResultsToFile
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`bool VtkOutput::IsCompatibleVariable(const std::string& rVariableName) const`

# Called by

- [WriteNodalResultsToFile](../../../../functions/input_output/vtk_output/VtkOutput/WriteNodalResultsToFile.md)
- [WriteElementResultsToFile](../../../../functions/input_output/vtk_output/VtkOutput/WriteElementResultsToFile.md)
- [WriteConditionResultsToFile](../../../../functions/input_output/vtk_output/VtkOutput/WriteConditionResultsToFile.md)