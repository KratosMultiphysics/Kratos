---
type: C++ Method
title: IdentifyFaces
resource: mpi/utilities/mpi_normal_calculation_utilities.cpp#L222-L256
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/mpi/utilities/mpi_normal_calculation_utilities/MPINormalCalculationUtils/CalculateOnSimplex
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void MPINormalCalculationUtils::IdentifyFaces(ModelPart& rModelPart, const Variable<double>& rVariable, int& MaxNeigh, int& NumNodes)`

# Calls

- [GetGeometry](../../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)

# Called by

- [CalculateOnSimplex](../../../../../functions/mpi/utilities/mpi_normal_calculation_utilities/MPINormalCalculationUtils/CalculateOnSimplex.md)