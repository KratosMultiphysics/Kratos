---
type: C++ Module
title: time_averaging_process
resource: processes/time_averaging_process.cpp#L1-L470
generated:
  by: okf-rs/0.7.0
relationships:
  imports:
  - target: external/sstream
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-define-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-variable-utils-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/time-averaging-process-h
    resolved_by: tree-sitter
    confidence: exact
  member_of:
  - target: packages/KratosMultiphysics
    resolved_by: tree-sitter
    confidence: exact
---

# Contains

- [TimeAveragingProcess](../../functions/processes/time_averaging_process/TimeAveragingProcess/TimeAveragingProcess.md)
- [Check](../../functions/processes/time_averaging_process/TimeAveragingProcess/Check.md)
- [ExecuteInitialize](../../functions/processes/time_averaging_process/TimeAveragingProcess/ExecuteInitialize.md)
- [ExecuteFinalizeSolutionStep](../../functions/processes/time_averaging_process/TimeAveragingProcess/ExecuteFinalizeSolutionStep.md)
- [if](../../functions/processes/time_averaging_process/if.md)
- [GetDefaultParameters](../../functions/processes/time_averaging_process/TimeAveragingProcess/GetDefaultParameters.md)
- [IsIntegrationStep](../../functions/processes/time_averaging_process/TimeAveragingProcess/IsIntegrationStep.md)
- [Info](../../functions/processes/time_averaging_process/TimeAveragingProcess/Info.md)
- [PrintInfo](../../functions/processes/time_averaging_process/TimeAveragingProcess/PrintInfo.md)
- [PrintData](../../functions/processes/time_averaging_process/TimeAveragingProcess/PrintData.md)
- [GetTimeAveragingMethod](../../functions/processes/time_averaging_process/TimeAveragingProcess/GetTimeAveragingMethod.md)
- [CalculateTimeIntegratedNodalHistoricalQuantity](../../functions/processes/time_averaging_process/TimeAveragingProcess/CalculateTimeIntegratedNodalHistoricalQuantity.md)
- [CalculateTimeIntegratedNonHistoricalQuantity](../../functions/processes/time_averaging_process/TimeAveragingProcess/CalculateTimeIntegratedNonHistoricalQuantity.md)
- [AverageMethod](../../functions/processes/time_averaging_process/TimeAveragingProcess/AverageMethod.md)
- [RootMeanSquareMethod](../../functions/processes/time_averaging_process/TimeAveragingProcess/RootMeanSquareMethod.md)

# Imports

- `<sstream>`
- `"includes/define.h"`
- `"utilities/variable_utils.h"`
- `"time_averaging_process.h"`

# Member of

- [KratosMultiphysics](../../packages/KratosMultiphysics.md)