---
type: C++ Method
title: SumForceAndTorque
resource: utilities/force_and_torque_utils.cpp#L42-L80
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/tests/test_force_and_torque_utils/TestForceAndTorqueUtils/test_with_moment
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`std::array<array_1d<double,3>,2> ForceAndTorqueUtils::SumForceAndTorque( const ModelPart& rModelPart, const Variable<array_1d<double,3>>& rForceVariable, const Variable<array_1d<double,3>>& rTorqueVariable)`

# Called by

- [test_with_moment](../../../../functions/tests/test_force_and_torque_utils/TestForceAndTorqueUtils/test_with_moment.md)