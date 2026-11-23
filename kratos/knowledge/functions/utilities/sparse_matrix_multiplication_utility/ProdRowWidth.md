---
type: C++ Function
title: ProdRowWidth
resource: utilities/sparse_matrix_multiplication_utility.h#L1068-L1124
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/utilities/sparse_matrix_multiplication_utility/MatrixMultiplicationRMerge
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static TIndex ProdRowWidth( const TIndex* AColumn, const TIndex* AColumnEnd, const TIndex* BPtr, const TIndex* BColumn, TIndex* Tmp1Column, TIndex* Tmp2Column, TIndex* Tmp3Column )`

# Called by

- [MatrixMultiplicationRMerge](../../../functions/utilities/sparse_matrix_multiplication_utility/MatrixMultiplicationRMerge.md)