---
type: C++ Method
title: PMultigridBuilderAndSolver
resource: solving_strategies/builder_and_solvers/p_multigrid/p_multigrid_builder_and_solver.cpp#L532-L538
visibility: private
generated:
  by: okf-rs/0.7.0
---

# Signature

`PMultigridBuilderAndSolver<TSparse,TDense>::PMultigridBuilderAndSolver(const typename LinearSolverType::Pointer& pSolver, Parameters Settings) : Interface(pSolver), mpImpl(new Impl(this))`