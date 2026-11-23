---
type: C++ Method
title: CalculateOnSimplex
resource: mpi/utilities/mpi_normal_calculation_utilities.cpp#L154-L197
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/mpi/utilities/mpi_normal_calculation_utilities/MPINormalCalculationUtils/IdentifyFaces
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/utilities/mpi_normal_calculation_utilities/MPINormalCalculationUtils/InitializeNormalData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/utilities/mpi_normal_calculation_utilities/MPINormalCalculationUtils/DetectEdges
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/utilities/mpi_normal_calculation_utilities/MPINormalCalculationUtils/UpdateNodeNormals
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void MPINormalCalculationUtils::CalculateOnSimplex(ModelPart& rModelPart, int Dimension, const Variable<double>& rVariable, const double rAlpha)`

# Calls

- [IdentifyFaces](../../../../../functions/mpi/utilities/mpi_normal_calculation_utilities/MPINormalCalculationUtils/IdentifyFaces.md)
- [InitializeNormalData](../../../../../functions/mpi/utilities/mpi_normal_calculation_utilities/MPINormalCalculationUtils/InitializeNormalData.md)
- [DetectEdges](../../../../../functions/mpi/utilities/mpi_normal_calculation_utilities/MPINormalCalculationUtils/DetectEdges.md)
- [UpdateNodeNormals](../../../../../functions/mpi/utilities/mpi_normal_calculation_utilities/MPINormalCalculationUtils/UpdateNodeNormals.md)