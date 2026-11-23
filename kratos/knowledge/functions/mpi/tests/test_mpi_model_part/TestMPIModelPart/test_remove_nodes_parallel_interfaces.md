---
type: Python Method
title: test_remove_nodes_parallel_interfaces
resource: mpi/tests/test_mpi_model_part.py#L13-L65
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python_scripts/testing/utilities/ReadDistributedModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/communicator/Communicator/MyPID
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part/ModelPart/RemoveNodesFromAllLevels
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/communicator/Communicator/GetNumberOfColors
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def test_remove_nodes_parallel_interfaces(self):`

# Calls

- [ReadDistributedModelPart](../../../../../functions/python_scripts/testing/utilities/ReadDistributedModelPart.md)
- [MyPID](../../../../../functions/sources/communicator/Communicator/MyPID.md)
- [RemoveNodesFromAllLevels](../../../../../functions/sources/model_part/ModelPart/RemoveNodesFromAllLevels.md)
- [GetNumberOfColors](../../../../../functions/sources/communicator/Communicator/GetNumberOfColors.md)