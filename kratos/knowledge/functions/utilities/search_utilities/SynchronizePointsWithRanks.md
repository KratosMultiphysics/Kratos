---
type: C++ Function
title: SynchronizePointsWithRanks
resource: utilities/search_utilities.h#L745-L913
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/amatrix_interface/Matrix/fill
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/utilities/search_utilities/SynchronizePointsWithBoundingBox
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static std::vector<int> SynchronizePointsWithRanks( TPointIteratorType itPointBegin, TPointIteratorType itPointEnd, std::vector<double>& rAllPointsCoordinates, std::vector<SignedIndexType>& rAllPointsLocalIds, std::vector<IndexType>& rAllPointsGlobalIds, const DataCommunicator& rDataCommunicator, const int NumberOfPoints, const int TotalNumberOfPoints, const bool IndexItIsJustPosition = false )`

# Calls

- [fill](../../../functions/includes/amatrix_interface/Matrix/fill.md)

# Called by

- [SynchronizePointsWithBoundingBox](../../../functions/utilities/search_utilities/SynchronizePointsWithBoundingBox.md)