---
type: C++ Function
title: CopyModelPartFlaggedNodalHistoricalVarToNonHistoricalVar
resource: utilities/variable_utils.h#L321-L354
generated:
  by: okf-rs/0.7.0
---

# Signature

`void CopyModelPartFlaggedNodalHistoricalVarToNonHistoricalVar( const Variable<TDataType>& rOriginVariable, const Variable<TDataType>& rDestinationVariable, const ModelPart& rOriginModelPart, ModelPart& rDestinationModelPart, const Flags& rFlag, const bool CheckValue = true, const unsigned int ReadBufferStep = 0)`