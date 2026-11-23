---
type: C++ Method
title: WriteNodalVectorValues
resource: input_output/vtk_output.cpp#L815-L828
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/vtk_output/VtkOutput/WriteVectorSolutionStepVariable
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/vtk_output/VtkOutput/WriteVectorContainerVariable
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/input_output/vtk_output/VtkOutput/WriteNodalContainerResults
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void VtkOutput::WriteNodalVectorValues( const ModelPart::NodesContainerType& rNodes, const TVarType& rVariable, const bool IsHistoricalValue, std::ofstream& rFileStream) const`

# Calls

- [WriteVectorSolutionStepVariable](../../../../functions/input_output/vtk_output/VtkOutput/WriteVectorSolutionStepVariable.md)
- [WriteVectorContainerVariable](../../../../functions/input_output/vtk_output/VtkOutput/WriteVectorContainerVariable.md)

# Called by

- [WriteNodalContainerResults](../../../../functions/input_output/vtk_output/VtkOutput/WriteNodalContainerResults.md)