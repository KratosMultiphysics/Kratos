---
type: C++ Function
title: SynchronousPointSynchronizationWithBoundingBox
resource: utilities/search_utilities.h#L344-L366
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/search_utilities/CalculateNumberOfPoints
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/search_utilities/SynchronizePointsWithBoundingBox
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static void SynchronousPointSynchronizationWithBoundingBox( TPointIteratorType itPointBegin, TPointIteratorType itPointEnd, DistributedSearchInformation& rSearchInfo, const TBoundingBoxType& rBoundingBox, const double ThresholdBoundingBox, const DataCommunicator& rDataCommunicator, const bool IndexItIsJustPosition = false )`

# Calls

- [CalculateNumberOfPoints](../../../functions/utilities/search_utilities/CalculateNumberOfPoints.md)
- [SynchronizePointsWithBoundingBox](../../../functions/utilities/search_utilities/SynchronizePointsWithBoundingBox.md)