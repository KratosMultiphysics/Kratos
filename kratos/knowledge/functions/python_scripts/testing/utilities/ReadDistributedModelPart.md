---
type: Python Function
title: ReadDistributedModelPart
resource: python_scripts/testing/utilities.py#L56-L82
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/mpi/tests/test_mpi_model_part/TestMPIModelPart/test_remove_nodes_parallel_interfaces
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/testing/utilities/ReadModelPart
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def ReadDistributedModelPart(mdpa_file_name, model_part, importer_settings=None):`

# Called by

- [test_remove_nodes_parallel_interfaces](../../../../functions/mpi/tests/test_mpi_model_part/TestMPIModelPart/test_remove_nodes_parallel_interfaces.md)
- [ReadModelPart](../../../../functions/python_scripts/testing/utilities/ReadModelPart.md)