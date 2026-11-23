---
type: C++ Function
title: FindPointOnMeshSimplified
resource: utilities/binbased_fast_point_locator.h#L278-L291
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/processes/apply_periodic_boundary_condition_process/ApplyPeriodicConditionProcess/ApplyConstraintsForPeriodicConditions
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python/add_geometrical_utilities_to_python/VectorizedFindHelper
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python/add_geometrical_utilities_to_python/AddGeometricalUtilitiesToPython
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/utilities/test_binbased_fast_point_locator/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/rve_periodicity_utility/RVEPeriodicityUtility/AuxiliaryAssignPeriodicity
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`bool FindPointOnMeshSimplified( const array_1d<double, 3 >& rCoordinates, Vector& rNShapeFunction, typename EntityType::Pointer& pEntity, const SizeType MaxNumberOfResults = 1000, const double Tolerance = 1.0e-5 )`

# Called by

- [ApplyConstraintsForPeriodicConditions](../../../functions/processes/apply_periodic_boundary_condition_process/ApplyPeriodicConditionProcess/ApplyConstraintsForPeriodicConditions.md)
- [VectorizedFindHelper](../../../functions/python/add_geometrical_utilities_to_python/VectorizedFindHelper.md)
- [AddGeometricalUtilitiesToPython](../../../functions/python/add_geometrical_utilities_to_python/AddGeometricalUtilitiesToPython.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../functions/tests/cpp_tests/utilities/test_binbased_fast_point_locator/KRATOS_TEST_CASE_IN_SUITE.md)
- [AuxiliaryAssignPeriodicity](../../../functions/utilities/rve_periodicity_utility/RVEPeriodicityUtility/AuxiliaryAssignPeriodicity.md)