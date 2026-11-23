---
type: C++ Method
title: CalculateNormal2D
resource: utilities/normal_calculation_utils.cpp#L428-L441
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/utilities/normal_calculation_utils/NormalCalculationUtils/AuxiliaryCalculateOnSimplex
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/normal_calculation_utils/CalculateOnSimplex
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/normal_calculation_utils/CalculateOnSimplexLowMemory
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void NormalCalculationUtils::CalculateNormal2D( Condition& rCondition, array_1d<double,3>& rAn, const Variable<array_1d<double,3>>& rNormalVariable )`

# Calls

- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)

# Called by

- [AuxiliaryCalculateOnSimplex](../../../../functions/utilities/normal_calculation_utils/NormalCalculationUtils/AuxiliaryCalculateOnSimplex.md)
- [CalculateOnSimplex](../../../../functions/utilities/normal_calculation_utils/CalculateOnSimplex.md)
- [CalculateOnSimplexLowMemory](../../../../functions/utilities/normal_calculation_utils/CalculateOnSimplexLowMemory.md)