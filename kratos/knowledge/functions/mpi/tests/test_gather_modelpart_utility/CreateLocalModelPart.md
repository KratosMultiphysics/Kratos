---
type: Python Function
title: CreateLocalModelPart
resource: mpi/tests/test_gather_modelpart_utility.py#L64-L79
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/model_part/ModelPart/CreateNewProperties
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/tests/test_gather_modelpart_utility/TestMesh/LocalNodes
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/tests/test_gather_modelpart_utility/TestMesh/LocalElements
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/tests/test_gather_modelpart_utility/TestMesh/PartitionIndex
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/mpi/tests/test_gather_modelpart_utility/TestGatherModelPartUtility/setUp
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def CreateLocalModelPart(model):`

# Calls

- [CreateNewProperties](../../../../functions/sources/model_part/ModelPart/CreateNewProperties.md)
- [LocalNodes](../../../../functions/mpi/tests/test_gather_modelpart_utility/TestMesh/LocalNodes.md)
- [LocalElements](../../../../functions/mpi/tests/test_gather_modelpart_utility/TestMesh/LocalElements.md)
- [PartitionIndex](../../../../functions/mpi/tests/test_gather_modelpart_utility/TestMesh/PartitionIndex.md)

# Called by

- [setUp](../../../../functions/mpi/tests/test_gather_modelpart_utility/TestGatherModelPartUtility/setUp.md)