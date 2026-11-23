---
type: C++ Function
title: CopyModelPartFlaggedNodalNonHistoricalVarToHistoricalVar
resource: utilities/variable_utils.h#L434-L445
generated:
  by: okf-rs/0.7.0
---

# Signature

`void CopyModelPartFlaggedNodalNonHistoricalVarToHistoricalVar( const Variable<TDataType>& rOriginVariable, const Variable<TDataType>& rDestinationVariable, ModelPart& rModelPart, const Flags& rFlag, const bool CheckValue = true, const unsigned int WriteBufferStep = 0)`