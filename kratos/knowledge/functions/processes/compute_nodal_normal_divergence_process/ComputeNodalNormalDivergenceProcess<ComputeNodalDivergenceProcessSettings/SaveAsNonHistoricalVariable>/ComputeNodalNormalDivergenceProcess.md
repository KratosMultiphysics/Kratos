---
type: C++ Method
title: ComputeNodalNormalDivergenceProcess
resource: processes/compute_nodal_normal_divergence_process.cpp#L164-L185
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/variable_utils/CheckVariableExists
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`ComputeNodalNormalDivergenceProcess<ComputeNodalDivergenceProcessSettings::SaveAsNonHistoricalVariable>::ComputeNodalNormalDivergenceProcess( ModelPart& rModelPart, const Variable<array_1d<double,3>>& rOriginVariable, const Variable<double>& rDivergenceVariable, const Variable<double>& rAreaVariable, const bool NormalizeDivergence, const bool NonHistoricalOriginVariable ) : mrModelPart(rModelPart), mpOriginVariable(&rOriginVariable), mpDivergenceVariable(&rDivergenceVariable), mpAreaVariable(&rAreaVariable), mNormalizeDivergence(NormalizeDivergence), mNonHistoricalOriginVariable(NonHistoricalOriginVariable)`

# Calls

- [CheckVariableExists](../../../../../functions/utilities/variable_utils/CheckVariableExists.md)