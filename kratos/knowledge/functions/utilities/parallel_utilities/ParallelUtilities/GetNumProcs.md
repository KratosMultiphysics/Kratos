---
type: C++ Method
title: GetNumProcs
resource: utilities/parallel_utilities.cpp#L69-L85
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/utilities/openmp_utils/GetNumberOfProcessors
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/parallel_utilities/ParallelUtilities/SetNumThreads
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`int ParallelUtilities::GetNumProcs()`

# Called by

- [GetNumberOfProcessors](../../../../functions/utilities/openmp_utils/GetNumberOfProcessors.md)
- [SetNumThreads](../../../../functions/utilities/parallel_utilities/ParallelUtilities/SetNumThreads.md)