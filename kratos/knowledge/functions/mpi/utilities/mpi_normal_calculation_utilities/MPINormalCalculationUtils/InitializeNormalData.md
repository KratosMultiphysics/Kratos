---
type: C++ Method
title: InitializeNormalData
resource: mpi/utilities/mpi_normal_calculation_utilities.cpp#L262-L403
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/utilities/mpi_normal_calculation_utilities/MPINormalCalculationUtils/FaceNormal3D
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/communicator/Communicator/GetNumberOfColors
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/mpi/utilities/mpi_normal_calculation_utilities/MPINormalCalculationUtils/CalculateOnSimplex
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void MPINormalCalculationUtils::InitializeNormalData(ModelPart& rModelPart, const Variable<double>& rVariable, int MaxNeigh, std::vector<double>& rNormals, std::vector<int>& rActiveNeigh, const int NumNodes)`

# Calls

- [GetGeometry](../../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [FaceNormal3D](../../../../../functions/mpi/utilities/mpi_normal_calculation_utilities/MPINormalCalculationUtils/FaceNormal3D.md)
- [GetNumberOfColors](../../../../../functions/sources/communicator/Communicator/GetNumberOfColors.md)

# Called by

- [CalculateOnSimplex](../../../../../functions/mpi/utilities/mpi_normal_calculation_utilities/MPINormalCalculationUtils/CalculateOnSimplex.md)