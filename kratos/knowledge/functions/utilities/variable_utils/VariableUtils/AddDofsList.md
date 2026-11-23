---
type: C++ Method
title: AddDofsList
resource: utilities/variable_utils.cpp#L119-L141
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/mpi/tests/test_mpi_communicator/TestMPICommunicator/test_GlobalNumberOf_Methods
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_dofs/TestDofs/testDofArrayUtilities
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_variable_utils/TestVariableUtils/testAddDofsList
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_variable_utils/TestVariableUtils/testAddDofsWithReactionList
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void VariableUtils::AddDofsList( const std::vector<std::string>& rDofsVarNamesList, ModelPart& rModelPart)`

# Called by

- [test_GlobalNumberOf_Methods](../../../../functions/mpi/tests/test_mpi_communicator/TestMPICommunicator/test_GlobalNumberOf_Methods.md)
- [testDofArrayUtilities](../../../../functions/tests/test_dofs/TestDofs/testDofArrayUtilities.md)
- [testAddDofsList](../../../../functions/tests/test_variable_utils/TestVariableUtils/testAddDofsList.md)
- [testAddDofsWithReactionList](../../../../functions/tests/test_variable_utils/TestVariableUtils/testAddDofsWithReactionList.md)