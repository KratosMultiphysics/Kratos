---
type: C++ Method
title: NormalContribution
resource: mpi/utilities/mpi_normal_calculation_utilities.cpp#L526-L545
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/geometries/geometry/PointsNumber
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/mpi/utilities/mpi_normal_calculation_utilities/MPINormalCalculationUtils/OrientFaces
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void MPINormalCalculationUtils::NormalContribution(Geometry<Node > &rGeom)`

# Calls

- [PointsNumber](../../../../../functions/geometries/geometry/PointsNumber.md)

# Called by

- [OrientFaces](../../../../../functions/mpi/utilities/mpi_normal_calculation_utilities/MPINormalCalculationUtils/OrientFaces.md)