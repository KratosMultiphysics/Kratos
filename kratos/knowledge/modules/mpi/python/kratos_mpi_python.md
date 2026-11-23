---
type: C++ Module
title: kratos_mpi_python
resource: mpi/python/kratos_mpi_python.cpp#L1-L78
generated:
  by: okf-rs/0.7.0
relationships:
  imports:
  - target: external/pybind11-pybind11-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/mpi-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/mpi-includes-mpi-communicator-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/mpi-includes-mpi-data-communicator-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/mpi-utilities-parallel-fill-communicator-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/add-mpi-communicator-to-python-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/add-mpi-data-communicator-to-python-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/add-mpi-utilities-to-python-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/add-mpi-debug-utilities-to-python-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/add-distributed-sparse-matrices-to-python-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-parallel-environment-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/add-mpi-search-strategies-to-python-h
    resolved_by: tree-sitter
    confidence: exact
  member_of:
  - target: packages/KratosMultiphysics
    resolved_by: tree-sitter
    confidence: exact
---

# Contains

- [InitializeMPIParallelRun](../../../functions/mpi/python/kratos_mpi_python/InitializeMPIParallelRun.md)
- [PYBIND11_MODULE](../../../functions/mpi/python/kratos_mpi_python/PYBIND11_MODULE.md)

# Imports

- `<pybind11/pybind11.h>`
- `"mpi.h"`
- `"mpi/includes/mpi_communicator.h"`
- `"mpi/includes/mpi_data_communicator.h"`
- `"mpi/utilities/parallel_fill_communicator.h"`
- `"add_mpi_communicator_to_python.h"`
- `"add_mpi_data_communicator_to_python.h"`
- `"add_mpi_utilities_to_python.h"`
- `"add_mpi_debug_utilities_to_python.h"`
- `"add_distributed_sparse_matrices_to_python.h"`
- `"includes/parallel_environment.h"`
- `"add_mpi_search_strategies_to_python.h"`

# Member of

- [KratosMultiphysics](../../../packages/KratosMultiphysics.md)