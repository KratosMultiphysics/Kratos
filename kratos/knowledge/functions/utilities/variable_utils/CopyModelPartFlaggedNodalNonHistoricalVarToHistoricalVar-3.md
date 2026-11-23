---
type: C++ Function
title: CopyModelPartFlaggedNodalNonHistoricalVarToHistoricalVar
resource: utilities/variable_utils.h#L448-L459
generated:
  by: okf-rs/0.7.0
---

# Signature

`void CopyModelPartFlaggedNodalNonHistoricalVarToHistoricalVar( const Variable<TDataType>& rVariable, const ModelPart& rOriginModelPart, ModelPart& rDestinationModelPart, const Flags& rFlag, const bool CheckValue = true, const unsigned int WriteBufferStep = 0)`