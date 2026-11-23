---
type: C++ Method
title: UpdateNodeNormals
resource: mpi/utilities/mpi_normal_calculation_utilities.cpp#L579-L628
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/utilities/mpi_normal_calculation_utilities/MPINormalCalculationUtils/FaceNormal2D
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/utilities/mpi_normal_calculation_utilities/MPINormalCalculationUtils/FaceNormal3D
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/mpi/utilities/mpi_normal_calculation_utilities/MPINormalCalculationUtils/CalculateOnSimplex
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void MPINormalCalculationUtils::UpdateNodeNormals(ModelPart &rModelPart, const unsigned int Dimension, const Variable<double> &rVariable)`

# Calls

- [GetGeometry](../../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [FaceNormal2D](../../../../../functions/mpi/utilities/mpi_normal_calculation_utilities/MPINormalCalculationUtils/FaceNormal2D.md)
- [FaceNormal3D](../../../../../functions/mpi/utilities/mpi_normal_calculation_utilities/MPINormalCalculationUtils/FaceNormal3D.md)

# Called by

- [CalculateOnSimplex](../../../../../functions/mpi/utilities/mpi_normal_calculation_utilities/MPINormalCalculationUtils/CalculateOnSimplex.md)