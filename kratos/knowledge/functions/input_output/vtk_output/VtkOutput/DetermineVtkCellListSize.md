---
type: C++ Method
title: DetermineVtkCellListSize
resource: input_output/vtk_output.cpp#L395-L404
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
  - target: functions/geometries/geometry/PointsNumber
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/input_output/vtk_output/VtkOutput/WriteConditionsAndElementsToFile
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`std::size_t VtkOutput::DetermineVtkCellListSize(const TContainerType& rContainer) const`

# Calls

- [SkippableEntity](../../../../functions/utilities/input_output_utilities/InputOutputUtilities/SkippableEntity.md)
- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [PointsNumber](../../../../functions/geometries/geometry/PointsNumber.md)

# Called by

- [WriteConditionsAndElementsToFile](../../../../functions/input_output/vtk_output/VtkOutput/WriteConditionsAndElementsToFile.md)