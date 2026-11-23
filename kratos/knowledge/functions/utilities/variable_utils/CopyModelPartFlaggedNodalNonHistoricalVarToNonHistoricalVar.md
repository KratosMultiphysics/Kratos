---
type: C++ Function
title: CopyModelPartFlaggedNodalNonHistoricalVarToNonHistoricalVar
resource: utilities/variable_utils.h#L475-L504
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/model_part/FullName
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void CopyModelPartFlaggedNodalNonHistoricalVarToNonHistoricalVar( const Variable<TDataType>& rOriginVariable, const Variable<TDataType>& rDestinationVariable, const ModelPart& rOriginModelPart, ModelPart& rDestinationModelPart, const Flags& rFlag, const bool CheckValue = true)`

# Calls

- [FullName](../../../functions/includes/model_part/FullName.md)