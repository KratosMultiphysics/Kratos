---
type: C++ Method
title: SetUpEffectiveDofArray
resource: utilities/dof_utilities/dof_array_utilities.cpp#L102-L170
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/communicator/Communicator/MyPID
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/future/solving_strategies/schemes/implicit_scheme/ImplicitScheme/SetUpDofArrays
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/utilities/test_dof_array_utilities/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_dofs/TestDofs/testDofArrayUtilities
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void DofArrayUtilities::SetUpEffectiveDofArray( const ModelPart &rModelPart, const DofsArrayType &rDofArray, DofsArrayType &rEffectiveDofArray, const unsigned int EchoLevel)`

# Calls

- [MyPID](../../../../../functions/sources/communicator/Communicator/MyPID.md)

# Called by

- [SetUpDofArrays](../../../../../functions/future/solving_strategies/schemes/implicit_scheme/ImplicitScheme/SetUpDofArrays.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../../functions/tests/cpp_tests/utilities/test_dof_array_utilities/KRATOS_TEST_CASE_IN_SUITE.md)
- [testDofArrayUtilities](../../../../../functions/tests/test_dofs/TestDofs/testDofArrayUtilities.md)