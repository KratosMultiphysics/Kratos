---
type: C++ Function
title: FastProjectDirection
resource: utilities/geometrical_projection_utilities.h#L87-L117
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/amatrix_interface/norm_2
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/processes/simple_mortar_mapper_process/SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>/AssemblyMortarOperators
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/utilities/test_geometrical_projection_utilities/TestFastProjectDirection
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/exact_mortar_segmentation_utility/ExactMortarIntegrationUtility<2, 2, false>/GetExactIntegration
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/exact_mortar_segmentation_utility/ExactMortarIntegrationUtility<2, 2, true>/GetExactIntegration
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/geometrical_projection_utilities/FastProjectOnGeometry
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static inline double FastProjectDirection( const TGeometryType& rGeom, const Point& rPointToProject, Point& rPointProjected, const array_1d<double,3>& rNormal, const array_1d<double,3>& rVector, const SizeType EchoLevel = 0 )`

# Calls

- [norm_2](../../../functions/includes/amatrix_interface/norm_2.md)

# Called by

- [AssemblyMortarOperators](../../../functions/processes/simple_mortar_mapper_process/SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>/AssemblyMortarOperators.md)
- [TestFastProjectDirection](../../../functions/tests/cpp_tests/utilities/test_geometrical_projection_utilities/TestFastProjectDirection.md)
- [GetExactIntegration](../../../functions/utilities/exact_mortar_segmentation_utility/ExactMortarIntegrationUtility<2, 2, false>/GetExactIntegration.md)
- [GetExactIntegration](../../../functions/utilities/exact_mortar_segmentation_utility/ExactMortarIntegrationUtility<2, 2, true>/GetExactIntegration.md)
- [FastProjectOnGeometry](../../../functions/utilities/geometrical_projection_utilities/FastProjectOnGeometry.md)