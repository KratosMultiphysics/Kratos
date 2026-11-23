---
type: C++ Method
title: perform_matrix_scaling_row
resource: linear_solvers/scaling_solver.h#L373-L392
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/linear_solvers/scaling_solver/ScalingSolver/SymmetricScaling
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static void perform_matrix_scaling_row( std::size_t CurrentRowIndex, int number_of_entries_in_row, IndexIteratorType row_it, // Should be index1_data().begin() + CurrentRowIndex IndexIteratorType index2_begin, ValueIteratorType value_begin, const VectorType& weights )`

# Called by

- [SymmetricScaling](../../../../functions/linear_solvers/scaling_solver/ScalingSolver/SymmetricScaling.md)