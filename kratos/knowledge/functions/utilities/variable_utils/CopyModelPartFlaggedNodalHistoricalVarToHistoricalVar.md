---
type: C++ Function
title: CopyModelPartFlaggedNodalHistoricalVarToHistoricalVar
resource: utilities/variable_utils.h#L229-L280
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/model_part/FullName
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void CopyModelPartFlaggedNodalHistoricalVarToHistoricalVar( const Variable<TDataType>& rOriginVariable, const Variable<TDataType>& rDestinationVariable, const ModelPart& rOriginModelPart, ModelPart& rDestinationModelPart, const Flags& rFlag, const bool CheckValue = true, const unsigned int ReadBufferStep = 0, const unsigned int WriteBufferStep = 0)`

# Calls

- [FullName](../../../functions/includes/model_part/FullName.md)