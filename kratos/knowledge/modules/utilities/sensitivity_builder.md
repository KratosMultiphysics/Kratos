---
type: C++ Module
title: sensitivity_builder
resource: utilities/sensitivity_builder.cpp#L1-L638
generated:
  by: okf-rs/0.7.0
relationships:
  imports:
  - target: external/containers-global-pointers-vector-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-kratos-parameters-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-parallel-environment-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-ublas-interface-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/input-output-logger-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/solving-strategies-schemes-sensitivity-builder-scheme-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-openmp-utils-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-parallel-utilities-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-pointer-map-communicator-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-variable-utils-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-sensitivity-builder-h
    resolved_by: tree-sitter
    confidence: exact
  member_of:
  - target: packages/KratosMultiphysics
    resolved_by: tree-sitter
    confidence: exact
---

# Contains

- [ExecuteFunctor](../../functions/utilities/sensitivity_builder/ExecuteFunctor.md)
- [ExecuteFunctorInContainer](../../functions/utilities/sensitivity_builder/ExecuteFunctorInContainer.md)
- [SetNonHistoricalValueToZeroFunctor](../../classes/utilities/sensitivity_builder/SetNonHistoricalValueToZeroFunctor.md)
- [SetHistoricalValueToZeroFunctor](../../classes/utilities/sensitivity_builder/SetHistoricalValueToZeroFunctor.md)
- [HasSensitivityContributions](../../functions/utilities/sensitivity_builder/HasSensitivityContributions.md)
- [GetDataFromVector](../../functions/utilities/sensitivity_builder/GetDataFromVector.md)
- [GetDataFromVector](../../functions/utilities/sensitivity_builder/GetDataFromVector-2.md)
- [AssembleContainerContributions](../../functions/utilities/sensitivity_builder/AssembleContainerContributions.md)
- [CalculateNodalSolutionStepSensitivityFunctor](../../classes/utilities/sensitivity_builder/CalculateNodalSolutionStepSensitivityFunctor.md)
- [CalculateNonHistoricalSensitivitiesFunctor](../../classes/utilities/sensitivity_builder/CalculateNonHistoricalSensitivitiesFunctor.md)
- [ReplaceDeprecatedNameIfExists](../../functions/utilities/sensitivity_builder/ReplaceDeprecatedNameIfExists.md)
- [SensitivityBuilder](../../functions/utilities/sensitivity_builder/SensitivityBuilder/SensitivityBuilder.md)
- [SensitivityBuilder](../../functions/utilities/sensitivity_builder/SensitivityBuilder/SensitivityBuilder-2.md)
- [SetResponseFunction](../../functions/utilities/sensitivity_builder/SensitivityBuilder/SetResponseFunction.md)
- [Initialize](../../functions/utilities/sensitivity_builder/SensitivityBuilder/Initialize.md)
- [InitializeSolutionStep](../../functions/utilities/sensitivity_builder/SensitivityBuilder/InitializeSolutionStep.md)
- [UpdateSensitivities](../../functions/utilities/sensitivity_builder/SensitivityBuilder/UpdateSensitivities.md)
- [FinalizeSolutionStep](../../functions/utilities/sensitivity_builder/SensitivityBuilder/FinalizeSolutionStep.md)
- [Finalize](../../functions/utilities/sensitivity_builder/SensitivityBuilder/Finalize.md)
- [Clear](../../functions/utilities/sensitivity_builder/SensitivityBuilder/Clear.md)
- [ClearFlags](../../functions/utilities/sensitivity_builder/SensitivityBuilder/ClearFlags.md)
- [ClearSensitivities](../../functions/utilities/sensitivity_builder/SensitivityBuilder/ClearSensitivities.md)
- [CalculateNodalSolutionStepSensitivities](../../functions/utilities/sensitivity_builder/SensitivityBuilder/CalculateNodalSolutionStepSensitivities.md)
- [CalculateNodalSolutionStepSensitivities](../../functions/utilities/sensitivity_builder/SensitivityBuilder/CalculateNodalSolutionStepSensitivities-2.md)
- [CalculateNonHistoricalSensitivities](../../functions/utilities/sensitivity_builder/SensitivityBuilder/CalculateNonHistoricalSensitivities.md)
- [CalculateNonHistoricalSensitivities](../../functions/utilities/sensitivity_builder/SensitivityBuilder/CalculateNonHistoricalSensitivities-2.md)
- [CalculateNonHistoricalSensitivities](../../functions/utilities/sensitivity_builder/SensitivityBuilder/CalculateNonHistoricalSensitivities-3.md)
- [CalculateNonHistoricalSensitivities](../../functions/utilities/sensitivity_builder/SensitivityBuilder/CalculateNonHistoricalSensitivities-4.md)
- [GetVariableLists](../../functions/utilities/sensitivity_builder/SensitivityBuilder/GetVariableLists.md)

# Imports

- `"containers/global_pointers_vector.h"`
- `"includes/kratos_parameters.h"`
- `"includes/parallel_environment.h"`
- `"includes/ublas_interface.h"`
- `"input_output/logger.h"`
- `"solving_strategies/schemes/sensitivity_builder_scheme.h"`
- `"utilities/openmp_utils.h"`
- `"utilities/parallel_utilities.h"`
- `"utilities/pointer_map_communicator.h"`
- `"utilities/variable_utils.h"`
- `"utilities/sensitivity_builder.h"`

# Member of

- [KratosMultiphysics](../../packages/KratosMultiphysics.md)