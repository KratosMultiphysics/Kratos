---
type: C++ Method
title: FaceNormal2D
resource: mpi/utilities/mpi_normal_calculation_utilities.cpp#L549-L556
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/mpi/utilities/mpi_normal_calculation_utilities/MPINormalCalculationUtils/OrientFaces
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/utilities/mpi_normal_calculation_utilities/MPINormalCalculationUtils/UpdateNodeNormals
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void MPINormalCalculationUtils::FaceNormal2D(array_1d<double,3> &An, Geometry<Node > &rGeometry)`

# Called by

- [OrientFaces](../../../../../functions/mpi/utilities/mpi_normal_calculation_utilities/MPINormalCalculationUtils/OrientFaces.md)
- [UpdateNodeNormals](../../../../../functions/mpi/utilities/mpi_normal_calculation_utilities/MPINormalCalculationUtils/UpdateNodeNormals.md)