---
type: C++ Method
title: WriteConnectivity
resource: input_output/vtk_output.cpp#L410-L434
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/input_output_utilities/InputOutputUtilities/SkippableEntity
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/geometry/PointsNumber
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/vtk_output/VtkOutput/ReorderConnectivity
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

`void VtkOutput::WriteConnectivity(const TContainerType& rContainer, std::ofstream& rFileStream) const`

# Calls

- [SkippableEntity](../../../../functions/utilities/input_output_utilities/InputOutputUtilities/SkippableEntity.md)
- [PointsNumber](../../../../functions/geometries/geometry/PointsNumber.md)
- [ReorderConnectivity](../../../../functions/input_output/vtk_output/VtkOutput/ReorderConnectivity.md)
- [WriteScalarDataToFile](../../../../functions/input_output/vtk_output/WriteScalarDataToFile.md)

# Called by

- [WriteConditionsAndElementsToFile](../../../../functions/input_output/vtk_output/VtkOutput/WriteConditionsAndElementsToFile.md)