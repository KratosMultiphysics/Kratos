---
type: C++ Function
title: CopyModelPartFlaggedNodalHistoricalVarToNonHistoricalVar
resource: utilities/variable_utils.h#L371-L382
generated:
  by: okf-rs/0.7.0
---

# Signature

`void CopyModelPartFlaggedNodalHistoricalVarToNonHistoricalVar( const Variable<TDataType>& rVariable, const ModelPart& rOriginModelPart, ModelPart& rDestinationModelPart, const Flags& rFlag, const bool CheckValue = true, const unsigned int ReadBufferStep = 0)`