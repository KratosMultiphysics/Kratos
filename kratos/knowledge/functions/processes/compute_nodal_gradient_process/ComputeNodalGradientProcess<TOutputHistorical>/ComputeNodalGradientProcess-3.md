---
type: C++ Method
title: ComputeNodalGradientProcess
resource: processes/compute_nodal_gradient_process.cpp#L61-L91
visibility: private
generated:
  by: okf-rs/0.7.0
---

# Signature

`ComputeNodalGradientProcess<TOutputHistorical>::ComputeNodalGradientProcess( ModelPart& rModelPart, const Variable<double>& rOriginVariable, const Variable<array_1d<double,3> >& rGradientVariable, const Variable<double>& rAreaVariable, const bool NonHistoricalVariable ) : mrModelPart(rModelPart), mpOriginVariable(&rOriginVariable), mpGradientVariable(&rGradientVariable), mpAreaVariable(&rAreaVariable), mNonHistoricalVariable(NonHistoricalVariable)`