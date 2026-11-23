---
type: C++ Function
title: SingularValueDecomposition
resource: utilities/svd_utils.h#L101-L122
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/SetDouble
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/RecursivelyValidateAndAssignDefaults
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetDouble
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetInt
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/svd_utils/JacobiSingularValueDecomposition
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static inline std::size_t SingularValueDecomposition( const TInputMatrix& rInputMatrix, TUmatrix& rUMatrix, TSMatrix& rSMatrix, TVMatrix& rVMatrix, Parameters ThisParameters )`

# Calls

- [SetDouble](../../../functions/sources/kratos_parameters/Parameters/SetDouble.md)
- [RecursivelyValidateAndAssignDefaults](../../../functions/sources/kratos_parameters/Parameters/RecursivelyValidateAndAssignDefaults.md)
- [GetString](../../../functions/sources/kratos_parameters/Parameters/GetString.md)
- [GetDouble](../../../functions/sources/kratos_parameters/Parameters/GetDouble.md)
- [GetInt](../../../functions/sources/kratos_parameters/Parameters/GetInt.md)
- [JacobiSingularValueDecomposition](../../../functions/utilities/svd_utils/JacobiSingularValueDecomposition.md)