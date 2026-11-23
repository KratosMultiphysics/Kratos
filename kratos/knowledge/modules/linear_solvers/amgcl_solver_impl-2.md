---
type: C++ Module
title: amgcl_solver_impl
resource: linear_solvers/amgcl_solver_impl.hpp#L1-L786
generated:
  by: okf-rs/0.7.0
relationships:
  imports:
  - target: external/amgcl-adapter-crs-tuple-hpp
    resolved_by: tree-sitter
    confidence: exact
  - target: external/amgcl-adapter-ublas-hpp
    resolved_by: tree-sitter
    confidence: exact
  - target: external/amgcl-adapter-zero-copy-hpp
    resolved_by: tree-sitter
    confidence: exact
  - target: external/amgcl-adapter-block-matrix-hpp
    resolved_by: tree-sitter
    confidence: exact
  - target: external/amgcl-backend-builtin-hpp
    resolved_by: tree-sitter
    confidence: exact
  - target: external/amgcl-value-type-static-matrix-hpp
    resolved_by: tree-sitter
    confidence: exact
  - target: external/amgcl-make-solver-hpp
    resolved_by: tree-sitter
    confidence: exact
  - target: external/amgcl-amg-hpp
    resolved_by: tree-sitter
    confidence: exact
  - target: external/amgcl-coarsening-runtime-hpp
    resolved_by: tree-sitter
    confidence: exact
  - target: external/amgcl-relaxation-runtime-hpp
    resolved_by: tree-sitter
    confidence: exact
  - target: external/amgcl-solver-runtime-hpp
    resolved_by: tree-sitter
    confidence: exact
  - target: external/amgcl-preconditioner-runtime-hpp
    resolved_by: tree-sitter
    confidence: exact
  - target: external/amgcl-coarsening-rigid-body-modes-hpp
    resolved_by: tree-sitter
    confidence: exact
  - target: external/amgcl-backend-vexcl-hpp
    resolved_by: tree-sitter
    confidence: exact
  - target: external/amgcl-backend-vexcl-static-matrix-hpp
    resolved_by: tree-sitter
    confidence: exact
  - target: external/amgcl-mpi-util-hpp
    resolved_by: tree-sitter
    confidence: exact
  - target: external/amgcl-mpi-make-solver-hpp
    resolved_by: tree-sitter
    confidence: exact
  - target: external/amgcl-mpi-preconditioner-hpp
    resolved_by: tree-sitter
    confidence: exact
  - target: external/amgcl-mpi-solver-runtime-hpp
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-define-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/input-output-logger-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/linear-solvers-amgcl-solver-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/variant
    resolved_by: tree-sitter
    confidence: exact
  member_of:
  - target: packages/KratosMultiphysics
    resolved_by: tree-sitter
    confidence: exact
---

# Contains

- [AMGCLAdaptor](../../classes/linear_solvers/amgcl_solver_impl/AMGCLAdaptor.md)
- [register_vexcl_static_matrix_type](../../functions/linear_solvers/amgcl_solver_impl/register_vexcl_static_matrix_type.md)
- [AMGCLSolver](../../functions/linear_solvers/amgcl_solver_impl/AMGCLSolver<TSparse,TDense>/AMGCLSolver.md)
- [AMGCLSolver](../../functions/linear_solvers/amgcl_solver_impl/AMGCLSolver<TSparse,TDense>/AMGCLSolver-2.md)
- [AMGCLSolver](../../functions/linear_solvers/amgcl_solver_impl/AMGCLSolver<TSparse,TDense>/AMGCLSolver-3.md)
- [AMGCLSolver](../../functions/linear_solvers/amgcl_solver_impl/AMGCLSolver<TSparse,TDense>/AMGCLSolver-4.md)
- [AMGCLSolver](../../functions/linear_solvers/amgcl_solver_impl/AMGCLSolver<TSparse,TDense>/AMGCLSolver-5.md)
- [ApplySettings](../../functions/linear_solvers/amgcl_solver_impl/AMGCLSolver<TSparse,TDense>/ApplySettings.md)
- [ProvideAdditionalData](../../functions/linear_solvers/amgcl_solver_impl/AMGCLSolver<TSparse,TDense>/ProvideAdditionalData.md)
- [Clear](../../functions/linear_solvers/amgcl_solver_impl/AMGCLSolver<TSparse,TDense>/Clear.md)
- [InitializeSolutionStep](../../functions/linear_solvers/amgcl_solver_impl/AMGCLSolver<TSparse,TDense>/InitializeSolutionStep.md)
- [AMGCLStaticVectorTraits](../../classes/linear_solvers/amgcl_solver_impl/AMGCLStaticVectorTraits.md)
- [PerformSolutionStep](../../functions/linear_solvers/amgcl_solver_impl/AMGCLSolver<TSparse,TDense>/PerformSolutionStep.md)

# Imports

- `<amgcl/adapter/crs_tuple.hpp>`
- `<amgcl/adapter/ublas.hpp>`
- `<amgcl/adapter/zero_copy.hpp>`
- `<amgcl/adapter/block_matrix.hpp>`
- `<amgcl/backend/builtin.hpp>`
- `<amgcl/value_type/static_matrix.hpp>`
- `<amgcl/make_solver.hpp>`
- `<amgcl/amg.hpp>`
- `<amgcl/coarsening/runtime.hpp>`
- `<amgcl/relaxation/runtime.hpp>`
- `<amgcl/solver/runtime.hpp>`
- `<amgcl/preconditioner/runtime.hpp>`
- `<amgcl/coarsening/rigid_body_modes.hpp>`
- `<amgcl/backend/vexcl.hpp>`
- `<amgcl/backend/vexcl_static_matrix.hpp>`
- `<amgcl/mpi/util.hpp>`
- `<amgcl/mpi/make_solver.hpp>`
- `<amgcl/mpi/preconditioner.hpp>`
- `<amgcl/mpi/solver/runtime.hpp>`
- `"includes/define.h"`
- `"input_output/logger.h"`
- `"linear_solvers/amgcl_solver.h"`
- `<variant>`

# Member of

- [KratosMultiphysics](../../packages/KratosMultiphysics.md)