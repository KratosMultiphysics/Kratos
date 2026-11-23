---
type: C++ Function
title: CopyModelPartFlaggedNodalNonHistoricalVarToHistoricalVar
resource: utilities/variable_utils.h#L398-L431
generated:
  by: okf-rs/0.7.0
---

# Signature

`void CopyModelPartFlaggedNodalNonHistoricalVarToHistoricalVar( const Variable<TDataType>& rOriginVariable, const Variable<TDataType>& rDestinationVariable, const ModelPart& rOriginModelPart, ModelPart& rDestinationModelPart, const Flags& rFlag, const bool CheckValue = true, const unsigned int WriteBufferStep = 0)`