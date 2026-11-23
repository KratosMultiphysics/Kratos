---
type: C++ Function
title: SynchronizePointsWithBoundingBox
resource: utilities/search_utilities.h#L932-L1022
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/search_utilities/DistributedSearchInformation/Reserve
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/search_utilities/SynchronizePointsWithRanks
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/search_utilities/DistributedSearchInformation/Shrink
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/utilities/search_utilities/SynchronousPointSynchronizationWithBoundingBox
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static void SynchronizePointsWithBoundingBox( TPointIteratorType itPointBegin, TPointIteratorType itPointEnd, DistributedSearchInformation& rSearchInfo, const TBoundingBoxType& rBoundingBox, const double ThresholdBoundingBox, const DataCommunicator& rDataCommunicator, const int NumberOfPoints, const int TotalNumberOfPoints, const bool IndexItIsJustPosition = false )`

# Calls

- [Reserve](../../../functions/utilities/search_utilities/DistributedSearchInformation/Reserve.md)
- [SynchronizePointsWithRanks](../../../functions/utilities/search_utilities/SynchronizePointsWithRanks.md)
- [Shrink](../../../functions/utilities/search_utilities/DistributedSearchInformation/Shrink.md)

# Called by

- [SynchronousPointSynchronizationWithBoundingBox](../../../functions/utilities/search_utilities/SynchronousPointSynchronizationWithBoundingBox.md)