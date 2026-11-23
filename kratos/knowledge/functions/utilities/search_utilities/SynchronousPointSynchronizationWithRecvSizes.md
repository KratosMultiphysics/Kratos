---
type: C++ Function
title: SynchronousPointSynchronizationWithRecvSizes
resource: utilities/search_utilities.h#L380-L397
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/search_utilities/CalculateNumberOfPoints
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static std::vector<int> SynchronousPointSynchronizationWithRecvSizes( TPointIteratorType itPointBegin, TPointIteratorType itPointEnd, std::vector<double>& rAllPointsCoordinates, std::vector<SignedIndexType>& rAllPointsIds, const DataCommunicator& rDataCommunicator )`

# Calls

- [CalculateNumberOfPoints](../../../functions/utilities/search_utilities/CalculateNumberOfPoints.md)