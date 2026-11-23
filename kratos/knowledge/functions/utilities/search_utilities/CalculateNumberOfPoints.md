---
type: C++ Function
title: CalculateNumberOfPoints
resource: utilities/search_utilities.h#L1073-L1103
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/utilities/search_utilities/SynchronousPointSynchronization
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/search_utilities/SynchronousPointSynchronizationWithBoundingBox
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/search_utilities/SynchronousPointSynchronizationWithRecvSizes
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/search_utilities/SynchronousPointSynchronizationWithRadius
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static void CalculateNumberOfPoints( TPointIteratorType itPointBegin, TPointIteratorType itPointEnd, int& rNumberOfPoints, int& rTotalNumberOfPoints, const DataCommunicator& rDataCommunicator )`

# Called by

- [SynchronousPointSynchronization](../../../functions/utilities/search_utilities/SynchronousPointSynchronization.md)
- [SynchronousPointSynchronizationWithBoundingBox](../../../functions/utilities/search_utilities/SynchronousPointSynchronizationWithBoundingBox.md)
- [SynchronousPointSynchronizationWithRecvSizes](../../../functions/utilities/search_utilities/SynchronousPointSynchronizationWithRecvSizes.md)
- [SynchronousPointSynchronizationWithRadius](../../../functions/utilities/search_utilities/SynchronousPointSynchronizationWithRadius.md)