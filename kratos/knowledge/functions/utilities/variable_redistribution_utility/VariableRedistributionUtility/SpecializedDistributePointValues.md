---
type: C++ Method
title: SpecializedDistributePointValues
resource: utilities/variable_redistribution_utility.cpp#L413-L456
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/variable_redistribution_utility/VariableRedistributionUtility/ComputeNodalSizes
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void VariableRedistributionUtility::SpecializedDistributePointValues( ModelPart& rModelPart, TContainerType& rEntitiesContainer, const Variable< TValueType >& rPointVariable, const Variable< TValueType >& rDistributedVariable, double Tolerance, unsigned int MaximumIterations)`

# Calls

- [ComputeNodalSizes](../../../../functions/utilities/variable_redistribution_utility/VariableRedistributionUtility/ComputeNodalSizes.md)