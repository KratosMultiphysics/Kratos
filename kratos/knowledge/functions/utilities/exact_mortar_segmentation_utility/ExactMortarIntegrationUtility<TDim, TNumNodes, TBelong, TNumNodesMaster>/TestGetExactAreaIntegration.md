---
type: C++ Method
title: TestGetExactAreaIntegration
resource: utilities/exact_mortar_segmentation_utility.cpp#L1104-L1116
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/exact_mortar_segmentation_utility/ExactMortarIntegrationUtility<TDim, TNumNodes, TBelong, TNumNodesMaster>/GetExactAreaIntegration
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`double ExactMortarIntegrationUtility<TDim, TNumNodes, TBelong, TNumNodesMaster>::TestGetExactAreaIntegration( Condition::Pointer pSlaveCond, Condition::Pointer pMasterCond )`

# Calls

- [GetExactAreaIntegration](../../../../functions/utilities/exact_mortar_segmentation_utility/ExactMortarIntegrationUtility<TDim, TNumNodes, TBelong, TNumNodesMaster>/GetExactAreaIntegration.md)
- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)