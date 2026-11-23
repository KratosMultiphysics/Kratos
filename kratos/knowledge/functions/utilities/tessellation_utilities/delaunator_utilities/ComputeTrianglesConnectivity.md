---
type: C++ Function
title: ComputeTrianglesConnectivity
resource: utilities/tessellation_utilities/delaunator_utilities.cpp#L144-L196
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/tessellation_utilities/delaunator_utilities/InitializeTriangulateIO
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/tessellation_utilities/delaunator_utilities/CleanTriangulateIO
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`std::vector<std::size_t> ComputeTrianglesConnectivity(const std::vector<double>& rCoordinates)`

# Calls

- [InitializeTriangulateIO](../../../../functions/utilities/tessellation_utilities/delaunator_utilities/InitializeTriangulateIO.md)
- [CleanTriangulateIO](../../../../functions/utilities/tessellation_utilities/delaunator_utilities/CleanTriangulateIO.md)