---
type: C++ Function
title: ProdRow
resource: utilities/sparse_matrix_multiplication_utility.h#L1148-L1257
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/utilities/sparse_matrix_multiplication_utility/MatrixMultiplicationRMerge
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static void ProdRow( const TAIndex* AColumn, const TAIndex* AColumnEnd, const TAValue *AValue, const TBIndex* BPtr, const TBIndex* BColumn, const TBValue *BValue, TOutIndex* OutColumn, TOutValue *OutValue, TOutIndex* Tmp2Column, TOutValue *Tmp2Value, TOutIndex* Tmp3Column, TOutValue *Tmp3Value )`

# Called by

- [MatrixMultiplicationRMerge](../../../functions/utilities/sparse_matrix_multiplication_utility/MatrixMultiplicationRMerge.md)