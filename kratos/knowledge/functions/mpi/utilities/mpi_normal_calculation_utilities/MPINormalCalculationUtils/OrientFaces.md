---
type: C++ Method
title: OrientFaces
resource: mpi/utilities/mpi_normal_calculation_utilities.cpp#L60-L150
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/utilities/mpi_normal_calculation_utilities/MPINormalCalculationUtils/OrientElement
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/utilities/mpi_normal_calculation_utilities/MPINormalCalculationUtils/NormalContribution
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/utilities/mpi_normal_calculation_utilities/MPINormalCalculationUtils/FaceNormal3D
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/utilities/mpi_normal_calculation_utilities/MPINormalCalculationUtils/FaceNormal2D
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/geometry/PointsNumber
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void MPINormalCalculationUtils::OrientFaces(ModelPart& rModelPart, bool OutwardsPositive)`

# Calls

- [GetGeometry](../../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [OrientElement](../../../../../functions/mpi/utilities/mpi_normal_calculation_utilities/MPINormalCalculationUtils/OrientElement.md)
- [NormalContribution](../../../../../functions/mpi/utilities/mpi_normal_calculation_utilities/MPINormalCalculationUtils/NormalContribution.md)
- [FaceNormal3D](../../../../../functions/mpi/utilities/mpi_normal_calculation_utilities/MPINormalCalculationUtils/FaceNormal3D.md)
- [FaceNormal2D](../../../../../functions/mpi/utilities/mpi_normal_calculation_utilities/MPINormalCalculationUtils/FaceNormal2D.md)
- [PointsNumber](../../../../../functions/geometries/geometry/PointsNumber.md)