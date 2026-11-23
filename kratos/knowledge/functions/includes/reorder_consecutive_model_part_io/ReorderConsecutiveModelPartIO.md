---
type: C++ Function
title: ReorderConsecutiveModelPartIO
resource: includes/reorder_consecutive_model_part_io.h#L88-L90
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/mpi/python_scripts/distributed_import_model_part_utility/DistributedImportModelPartUtility/ImportModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/single_import_model_part/SingleImportModelPart/Import
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`ReorderConsecutiveModelPartIO( Kratos::shared_ptr<std::iostream> Stream, const Flags Options) = delete;`

# Called by

- [ImportModelPart](../../../functions/mpi/python_scripts/distributed_import_model_part_utility/DistributedImportModelPartUtility/ImportModelPart.md)
- [Import](../../../functions/utilities/single_import_model_part/SingleImportModelPart/Import.md)