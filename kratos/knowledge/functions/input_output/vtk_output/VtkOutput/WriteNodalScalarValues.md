---
type: C++ Method
title: WriteNodalScalarValues
resource: input_output/vtk_output.cpp#L796-L809
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/vtk_output/VtkOutput/WriteScalarSolutionStepVariable
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/vtk_output/VtkOutput/WriteScalarContainerVariable
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/input_output/vtk_output/VtkOutput/WriteNodalContainerResults
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void VtkOutput::WriteNodalScalarValues( const ModelPart::NodesContainerType& rNodes, const TVarType& rVariable, const bool IsHistoricalValue, std::ofstream& rFileStream) const`

# Calls

- [WriteScalarSolutionStepVariable](../../../../functions/input_output/vtk_output/VtkOutput/WriteScalarSolutionStepVariable.md)
- [WriteScalarContainerVariable](../../../../functions/input_output/vtk_output/VtkOutput/WriteScalarContainerVariable.md)

# Called by

- [WriteNodalContainerResults](../../../../functions/input_output/vtk_output/VtkOutput/WriteNodalContainerResults.md)