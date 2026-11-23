---
type: C++ Function
title: ComputePredictComponent
resource: solving_strategies/schemes/residual_based_bdf_custom_scheme.h#L523-L557
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/node/HasDofFor
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/solving_strategies/schemes/residual_based_bdf_custom_scheme/ResidualBasedBDFCustomScheme/Predict
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ComputePredictComponent( NodesArrayType::iterator itNode, const TClassVar& rVariable, const TClassVar& rDerivedVariable, const TClassVar& rDerived2Variable, const double DeltaTime )`

# Calls

- [HasDofFor](../../../../functions/includes/node/HasDofFor.md)

# Called by

- [Predict](../../../../functions/solving_strategies/schemes/residual_based_bdf_custom_scheme/ResidualBasedBDFCustomScheme/Predict.md)