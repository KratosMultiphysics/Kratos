---
type: C++ Function
title: partial_product_no_add
resource: spaces/ublas_space.h#L1068-L1095
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/spaces/ublas_space/ParallelProductNoAdd
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static void partial_product_no_add( int number_of_rows, typename compressed_matrix<TDataType>::index_array_type::const_iterator row_begin, typename compressed_matrix<TDataType>::index_array_type::const_iterator index2_begin, typename compressed_matrix<TDataType>::value_array_type::const_iterator value_begin, const VectorType& input_vec, unsigned int output_begin_index, VectorType& output_vec // typename VectorType::iterator output_vec_begin )`

# Called by

- [ParallelProductNoAdd](../../../functions/spaces/ublas_space/ParallelProductNoAdd.md)