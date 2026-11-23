---
type: C++ Function
title: GetPybindArray
resource: python/numpy_utils.h#L179-L212
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/containers/nd_data/NDData<TDataType>/pData
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`pybind11::array_t<TDataType, pybind11::array::c_style> GetPybindArray(NDData<TDataType>& rDDArray)`

# Calls

- [pData](../../../functions/containers/nd_data/NDData<TDataType>/pData.md)