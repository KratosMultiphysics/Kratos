---
type: C++ Method
title: Execute
resource: processes/compute_nodal_gradient_process.cpp#L97-L113
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/processes/compute_nodal_gradient_process/ComputeNodalGradientProcess<TOutputHistorical>/ClearGradient
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/compute_nodal_gradient_process/ComputeNodalGradientProcess<TOutputHistorical>/ComputeElementalContributionsAndVolume
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/compute_nodal_gradient_process/ComputeNodalGradientProcess<TOutputHistorical>/SynchronizeGradientAndVolume
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/compute_nodal_gradient_process/ComputeNodalGradientProcess<TOutputHistorical>/PonderateGradient
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ComputeNodalGradientProcess<TOutputHistorical>::Execute()`

# Calls

- [ClearGradient](../../../../functions/processes/compute_nodal_gradient_process/ComputeNodalGradientProcess<TOutputHistorical>/ClearGradient.md)
- [ComputeElementalContributionsAndVolume](../../../../functions/processes/compute_nodal_gradient_process/ComputeNodalGradientProcess<TOutputHistorical>/ComputeElementalContributionsAndVolume.md)
- [SynchronizeGradientAndVolume](../../../../functions/processes/compute_nodal_gradient_process/ComputeNodalGradientProcess<TOutputHistorical>/SynchronizeGradientAndVolume.md)
- [PonderateGradient](../../../../functions/processes/compute_nodal_gradient_process/ComputeNodalGradientProcess<TOutputHistorical>/PonderateGradient.md)