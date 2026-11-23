---
type: C++ Method
title: ComputeEquivalentForceAndTorque
resource: utilities/force_and_torque_utils.cpp#L83-L154
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/tests/test_force_and_torque_utils/TestForceAndTorqueUtils/test_with_moment
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_force_and_torque_utils/TestForceAndTorqueUtils/test_without_moment
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`std::array<array_1d<double,3>,2> ForceAndTorqueUtils::ComputeEquivalentForceAndTorque( const ModelPart& rModelPart, const array_1d<double,3>& rReferencePoint, const Variable<array_1d<double,3>>& rForceVariable, const Variable<array_1d<double,3>>& rTorqueVariable)`

# Called by

- [test_with_moment](../../../../functions/tests/test_force_and_torque_utils/TestForceAndTorqueUtils/test_with_moment.md)
- [test_without_moment](../../../../functions/tests/test_force_and_torque_utils/TestForceAndTorqueUtils/test_without_moment.md)