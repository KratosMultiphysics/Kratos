---
type: C++ Module
title: add_tensor_adaptors_to_python
resource: python/add_tensor_adaptors_to_python.cpp#L1-L345
generated:
  by: okf-rs/0.7.0
relationships:
  imports:
  - target: external/cstdint
    resolved_by: tree-sitter
    confidence: exact
  - target: external/pybind11-functional-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/pybind11-numpy-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/pybind11-stl-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-define-python-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-model-part-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/python-numpy-utils-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-container-io-utils-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-parallel-utilities-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/tensor-adaptors-tensor-adaptor-utils-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/tensor-adaptors-combined-tensor-adaptor-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/tensor-adaptors-connectivity-ids-tensor-adaptor-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/tensor-adaptors-equation-ids-tensor-adaptor-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/tensor-adaptors-fixity-tensor-adaptor-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/tensor-adaptors-flags-tensor-adaptor-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/tensor-adaptors-gauss-point-variable-tensor-adaptor-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/tensor-adaptors-historical-variable-tensor-adaptor-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/tensor-adaptors-node-position-tensor-adaptor-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/tensor-adaptors-geometries-tensor-adaptor-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/tensor-adaptors-tensor-adaptor-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/tensor-adaptors-variable-tensor-adaptor-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/tensor-adaptors-geometry-metrics-tensor-adaptor-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/add-tensor-adaptors-to-python-h
    resolved_by: tree-sitter
    confidence: exact
  member_of:
  - target: packages/KratosMultiphysics
    resolved_by: tree-sitter
    confidence: exact
---

# Contains

- [AddBaseTensorAdaptor](../../functions/python/add_tensor_adaptors_to_python/AddBaseTensorAdaptor.md)
- [AddCombinedTensorAdaptor](../../functions/python/add_tensor_adaptors_to_python/AddCombinedTensorAdaptor.md)
- [AddTensorAdaptorsToPython](../../functions/python/add_tensor_adaptors_to_python/AddTensorAdaptorsToPython.md)

# Imports

- `<cstdint>`
- `<pybind11/functional.h>`
- `<pybind11/numpy.h>`
- `<pybind11/stl.h>`
- `"includes/define_python.h"`
- `"includes/model_part.h"`
- `"python/numpy_utils.h"`
- `"utilities/container_io_utils.h"`
- `"utilities/parallel_utilities.h"`
- `"tensor_adaptors/tensor_adaptor_utils.h"`
- `"tensor_adaptors/combined_tensor_adaptor.h"`
- `"tensor_adaptors/connectivity_ids_tensor_adaptor.h"`
- `"tensor_adaptors/equation_ids_tensor_adaptor.h"`
- `"tensor_adaptors/fixity_tensor_adaptor.h"`
- `"tensor_adaptors/flags_tensor_adaptor.h"`
- `"tensor_adaptors/gauss_point_variable_tensor_adaptor.h"`
- `"tensor_adaptors/historical_variable_tensor_adaptor.h"`
- `"tensor_adaptors/node_position_tensor_adaptor.h"`
- `"tensor_adaptors/geometries_tensor_adaptor.h"`
- `"tensor_adaptors/tensor_adaptor.h"`
- `"tensor_adaptors/variable_tensor_adaptor.h"`
- `"tensor_adaptors/geometry_metrics_tensor_adaptor.h"`
- `"add_tensor_adaptors_to_python.h"`

# Member of

- [KratosMultiphysics](../../packages/KratosMultiphysics.md)