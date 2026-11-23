---
type: C++ Method
title: AuxiliaryAssignPeriodicity
resource: utilities/rve_periodicity_utility.cpp#L48-L139
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/amatrix_interface/prod
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/binbased_fast_point_locator/FindPointOnMeshSimplified
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/amatrix_interface/norm_2
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void RVEPeriodicityUtility::AuxiliaryAssignPeriodicity( ModelPart& rMasterModelPart, ModelPart& rSlaveModelPart, const Matrix& rStrainTensor, const Vector& rDirection, const double SearchTolerance )`

# Calls

- [prod](../../../../functions/includes/amatrix_interface/prod.md)
- [FindPointOnMeshSimplified](../../../../functions/utilities/binbased_fast_point_locator/FindPointOnMeshSimplified.md)
- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [norm_2](../../../../functions/includes/amatrix_interface/norm_2.md)