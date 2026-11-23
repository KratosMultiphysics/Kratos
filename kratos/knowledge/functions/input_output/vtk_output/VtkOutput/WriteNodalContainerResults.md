---
type: C++ Method
title: WriteNodalContainerResults
resource: input_output/vtk_output.cpp#L681-L715
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/vtk_output/VtkOutput/WriteNodalScalarValues
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/vtk_output/VtkOutput/WriteNodalVectorValues
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/input_output/vtk_output/VtkOutput/WriteNodalResultsToFile
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void VtkOutput::WriteNodalContainerResults( const std::string& rVariableName, const ModelPart::NodesContainerType& rNodes, const bool IsHistoricalValue, std::ofstream& rFileStream) const`

# Calls

- [WriteNodalScalarValues](../../../../functions/input_output/vtk_output/VtkOutput/WriteNodalScalarValues.md)
- [WriteNodalVectorValues](../../../../functions/input_output/vtk_output/VtkOutput/WriteNodalVectorValues.md)

# Called by

- [WriteNodalResultsToFile](../../../../functions/input_output/vtk_output/VtkOutput/WriteNodalResultsToFile.md)