---
type: C++ Method
title: GetExactAreaIntegration
resource: utilities/exact_mortar_segmentation_utility.cpp#L1026-L1041
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/exact_mortar_segmentation_utility/ExactMortarIntegrationUtility<TDim, TNumNodes, TBelong, TNumNodesMaster>/GetTotalArea
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/cpp_tests/processes/test_mortar_mapper_process/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/exact_mortar_segmentation_utility/ExactMortarIntegrationUtility<TDim, TNumNodes, TBelong, TNumNodesMaster>/TestGetExactAreaIntegration
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`bool ExactMortarIntegrationUtility<TDim, TNumNodes, TBelong, TNumNodesMaster>::GetExactAreaIntegration( const GeometryType& rOriginalSlaveGeometry, const array_1d<double, 3>& rSlaveNormal, const GeometryType& rOriginalMasterGeometry, const array_1d<double, 3>& rMasterNormal, double& rArea )`

# Calls

- [GetTotalArea](../../../../functions/utilities/exact_mortar_segmentation_utility/ExactMortarIntegrationUtility<TDim, TNumNodes, TBelong, TNumNodesMaster>/GetTotalArea.md)

# Called by

- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/processes/test_mortar_mapper_process/KRATOS_TEST_CASE_IN_SUITE.md)
- [TestGetExactAreaIntegration](../../../../functions/utilities/exact_mortar_segmentation_utility/ExactMortarIntegrationUtility<TDim, TNumNodes, TBelong, TNumNodesMaster>/TestGetExactAreaIntegration.md)