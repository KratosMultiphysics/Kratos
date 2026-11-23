---
type: C++ Method
title: DetectEdges
resource: mpi/utilities/mpi_normal_calculation_utilities.cpp#L407-L503
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/mpi/utilities/mpi_normal_calculation_utilities/MPINormalCalculationUtils/CalculateOnSimplex
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void MPINormalCalculationUtils::DetectEdges(ModelPart &rModelPart, unsigned int Dimension, double MaxAngle, const std::vector<double>& rNormals, const std::vector<int>& rActiveNeigh, const int NumNodes, const int MaxNeighs)`

# Called by

- [CalculateOnSimplex](../../../../../functions/mpi/utilities/mpi_normal_calculation_utilities/MPINormalCalculationUtils/CalculateOnSimplex.md)