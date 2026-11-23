---
type: Python Method
title: __execute_test
resource: mpi/tests/test_distributed_import_model_part_utility.py#L27-L99
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/AddEmptyValue
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/SetString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/SetBool
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/communicator/Communicator/GlobalNumberOfNodes
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/communicator/Communicator/GlobalNumberOfElements
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/communicator/Communicator/GlobalNumberOfConditions
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/model_part/NumberOfSubModelParts
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part/ModelPart/HasSubModelPart
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/mpi/tests/test_distributed_import_model_part_utility/TestDistributedImportModelPartUtility/test_with_all_ranks
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/tests/test_distributed_import_model_part_utility/TestDistributedImportModelPartUtility/test_with_sub_comm
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/tests/test_distributed_import_model_part_utility/TestDistributedImportModelPartUtility/test_in_memory_with_all_ranks
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/tests/test_distributed_import_model_part_utility/TestDistributedImportModelPartUtility/test_in_memory_with_sub_comm
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def __execute_test(self, in_memory, all_ranks):`

# Calls

- [AddEmptyValue](../../../../../functions/sources/kratos_parameters/Parameters/AddEmptyValue.md)
- [SetString](../../../../../functions/sources/kratos_parameters/Parameters/SetString.md)
- [SetBool](../../../../../functions/sources/kratos_parameters/Parameters/SetBool.md)
- [GlobalNumberOfNodes](../../../../../functions/sources/communicator/Communicator/GlobalNumberOfNodes.md)
- [GlobalNumberOfElements](../../../../../functions/sources/communicator/Communicator/GlobalNumberOfElements.md)
- [GlobalNumberOfConditions](../../../../../functions/sources/communicator/Communicator/GlobalNumberOfConditions.md)
- [NumberOfSubModelParts](../../../../../functions/includes/model_part/NumberOfSubModelParts.md)
- [HasSubModelPart](../../../../../functions/sources/model_part/ModelPart/HasSubModelPart.md)

# Called by

- [test_with_all_ranks](../../../../../functions/mpi/tests/test_distributed_import_model_part_utility/TestDistributedImportModelPartUtility/test_with_all_ranks.md)
- [test_with_sub_comm](../../../../../functions/mpi/tests/test_distributed_import_model_part_utility/TestDistributedImportModelPartUtility/test_with_sub_comm.md)
- [test_in_memory_with_all_ranks](../../../../../functions/mpi/tests/test_distributed_import_model_part_utility/TestDistributedImportModelPartUtility/test_in_memory_with_all_ranks.md)
- [test_in_memory_with_sub_comm](../../../../../functions/mpi/tests/test_distributed_import_model_part_utility/TestDistributedImportModelPartUtility/test_in_memory_with_sub_comm.md)