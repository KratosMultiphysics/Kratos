---
type: C++ Method
title: WriteCellType
resource: input_output/vtk_output.cpp#L440-L459
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/input_output_utilities/InputOutputUtilities/SkippableEntity
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/vtk_output/WriteScalarDataToFile
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/input_output/vtk_output/VtkOutput/WriteConditionsAndElementsToFile
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void VtkOutput::WriteCellType(const TContainerType& rContainer, std::ofstream& rFileStream) const`

# Calls

- [SkippableEntity](../../../../functions/utilities/input_output_utilities/InputOutputUtilities/SkippableEntity.md)
- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [WriteScalarDataToFile](../../../../functions/input_output/vtk_output/WriteScalarDataToFile.md)

# Called by

- [WriteConditionsAndElementsToFile](../../../../functions/input_output/vtk_output/VtkOutput/WriteConditionsAndElementsToFile.md)