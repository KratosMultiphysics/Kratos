---
type: C++ Function
title: CopyModelPartFlaggedNodalHistoricalVarToHistoricalVar
resource: utilities/variable_utils.h#L283-L299
generated:
  by: okf-rs/0.7.0
---

# Signature

`void CopyModelPartFlaggedNodalHistoricalVarToHistoricalVar( const Variable<TDataType>& rOriginVariable, const Variable<TDataType>& rDestinationVariable, ModelPart& rModelPart, const Flags& rFlag, const bool CheckValue = true, const unsigned int ReadBufferStep = 0, const unsigned int WriteBufferStep = 0)`