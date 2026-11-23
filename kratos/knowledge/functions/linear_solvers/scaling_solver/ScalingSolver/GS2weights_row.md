---
type: C++ Method
title: GS2weights_row
resource: linear_solvers/scaling_solver.h#L413-L432
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/linear_solvers/scaling_solver/ScalingSolver/GetScalingWeights
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static void GS2weights_row( std::size_t CurrentRowIndex, int number_of_entries_in_row, ConstIndexIteratorType row_it, ConstIndexIteratorType index2_begin, // Not used in original GS2weights, but kept for consistency if needed later ConstValueIteratorType value_begin, VectorType& weights )`

# Called by

- [GetScalingWeights](../../../../functions/linear_solvers/scaling_solver/ScalingSolver/GetScalingWeights.md)