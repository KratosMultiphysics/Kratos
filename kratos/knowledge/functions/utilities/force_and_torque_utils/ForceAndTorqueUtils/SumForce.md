---
type: C++ Method
title: SumForce
resource: utilities/force_and_torque_utils.cpp#L28-L39
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/tests/test_force_and_torque_utils/TestForceAndTorqueUtils/test_force
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`array_1d<double,3> ForceAndTorqueUtils::SumForce( const ModelPart& rModelPart, const Variable<array_1d<double,3>>& rForceVariable)`

# Called by

- [test_force](../../../../functions/tests/test_force_and_torque_utils/TestForceAndTorqueUtils/test_force.md)