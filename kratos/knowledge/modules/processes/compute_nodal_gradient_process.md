---
type: C++ Module
title: compute_nodal_gradient_process
resource: processes/compute_nodal_gradient_process.cpp#L1-L415
generated:
  by: okf-rs/0.7.0
relationships:
  imports:
  - target: external/containers-model-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-variable-utils-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-geometry-utilities-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/processes-compute-nodal-gradient-process-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-parallel-utilities-h
    resolved_by: tree-sitter
    confidence: exact
  member_of:
  - target: packages/KratosMultiphysics
    resolved_by: tree-sitter
    confidence: exact
---

# Contains

- [ComputeNodalGradientProcess](../../functions/processes/compute_nodal_gradient_process/ComputeNodalGradientProcess<TOutputHistorical>/ComputeNodalGradientProcess.md)
- [ComputeNodalGradientProcess](../../functions/processes/compute_nodal_gradient_process/ComputeNodalGradientProcess<TOutputHistorical>/ComputeNodalGradientProcess-2.md)
- [ComputeNodalGradientProcess](../../functions/processes/compute_nodal_gradient_process/ComputeNodalGradientProcess<TOutputHistorical>/ComputeNodalGradientProcess-3.md)
- [Execute](../../functions/processes/compute_nodal_gradient_process/ComputeNodalGradientProcess<TOutputHistorical>/Execute.md)
- [ComputeElementalContributionsAndVolume](../../functions/processes/compute_nodal_gradient_process/ComputeNodalGradientProcess<TOutputHistorical>/ComputeElementalContributionsAndVolume.md)
- [TLSType](../../classes/processes/compute_nodal_gradient_process/TLSType.md)
- [PrepareMemberVariablesFromSettings](../../functions/processes/compute_nodal_gradient_process/ComputeNodalGradientProcess<TOutputHistorical>/PrepareMemberVariablesFromSettings.md)
- [CheckOriginAndAreaVariables](../../functions/processes/compute_nodal_gradient_process/ComputeNodalGradientProcess<TOutputHistorical>/CheckOriginAndAreaVariables.md)
- [ClearGradient](../../functions/processes/compute_nodal_gradient_process/ComputeNodalGradientProcess<TOutputHistorical>/ClearGradient.md)
- [PonderateGradient](../../functions/processes/compute_nodal_gradient_process/ComputeNodalGradientProcess<TOutputHistorical>/PonderateGradient.md)
- [SynchronizeGradientAndVolume](../../functions/processes/compute_nodal_gradient_process/ComputeNodalGradientProcess<TOutputHistorical>/SynchronizeGradientAndVolume.md)
- [GetDefaultParameters](../../functions/processes/compute_nodal_gradient_process/ComputeNodalGradientProcess<TOutputHistorical>/GetDefaultParameters.md)
- [GetVariableVector](../../functions/processes/compute_nodal_gradient_process/VariableVectorRetriever<TInputHistorical>/GetVariableVector.md)

# Imports

- `"containers/model.h"`
- `"utilities/variable_utils.h"`
- `"utilities/geometry_utilities.h"`
- `"processes/compute_nodal_gradient_process.h"`
- `"utilities/parallel_utilities.h"`

# Member of

- [KratosMultiphysics](../../packages/KratosMultiphysics.md)