---
type: C++ Function
title: CheckVariableExists
resource: utilities/variable_utils.h#L1066-L1079
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/processes/compute_nodal_gradient_process/ComputeNodalGradientProcess<TOutputHistorical>/ComputeNodalGradientProcess
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/compute_nodal_gradient_process/ComputeNodalGradientProcess<TOutputHistorical>/CheckOriginAndAreaVariables
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/compute_nodal_normal_divergence_process/ComputeNodalNormalDivergenceProcess<ComputeNodalDivergenceProcessSettings/SaveAsHistoricalVariable>/ComputeNodalNormalDivergenceProcess
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/compute_nodal_normal_divergence_process/ComputeNodalNormalDivergenceProcess<ComputeNodalDivergenceProcessSettings/SaveAsNonHistoricalVariable>/ComputeNodalNormalDivergenceProcess
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/variable_utils/ApplyFixity
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/variable_utils/ApplyVector
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/variable_utils/AddDofWithReaction
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`int CheckVariableExists( const TVarType& rVariable, const NodesContainerType& rNodes )`

# Called by

- [ComputeNodalGradientProcess](../../../functions/processes/compute_nodal_gradient_process/ComputeNodalGradientProcess<TOutputHistorical>/ComputeNodalGradientProcess.md)
- [CheckOriginAndAreaVariables](../../../functions/processes/compute_nodal_gradient_process/ComputeNodalGradientProcess<TOutputHistorical>/CheckOriginAndAreaVariables.md)
- [ComputeNodalNormalDivergenceProcess](../../../functions/processes/compute_nodal_normal_divergence_process/ComputeNodalNormalDivergenceProcess<ComputeNodalDivergenceProcessSettings/SaveAsHistoricalVariable>/ComputeNodalNormalDivergenceProcess.md)
- [ComputeNodalNormalDivergenceProcess](../../../functions/processes/compute_nodal_normal_divergence_process/ComputeNodalNormalDivergenceProcess<ComputeNodalDivergenceProcessSettings/SaveAsNonHistoricalVariable>/ComputeNodalNormalDivergenceProcess.md)
- [ApplyFixity](../../../functions/utilities/variable_utils/ApplyFixity.md)
- [ApplyVector](../../../functions/utilities/variable_utils/ApplyVector.md)
- [AddDofWithReaction](../../../functions/utilities/variable_utils/AddDofWithReaction.md)