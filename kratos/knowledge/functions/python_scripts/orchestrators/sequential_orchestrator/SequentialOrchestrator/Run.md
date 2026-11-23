---
type: Python Method
title: Run
resource: python_scripts/orchestrators/sequential_orchestrator.py#L24-L70
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python_scripts/orchestrators/sequential_orchestrator/SequentialOrchestrator/__GetStagesToCheckpointList
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/orchestrators/orchestrator/Orchestrator/GetProject
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/IsString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/orchestrators/sequential_orchestrator/SequentialOrchestrator/__GetExecutionList
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/orchestrators/orchestrator/Orchestrator/CheckStageSettings
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/orchestrators/orchestrator/Orchestrator/CreateStage
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/orchestrators/orchestrator/Orchestrator/RunCurrentStagePreprocess
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/project/Project/GetOutputData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/analysis_stage/AnalysisStage/GetFinalData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/orchestrators/orchestrator/Orchestrator/RunCurrentStagePostprocess
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/orchestrators/sequential_orchestrator/SequentialOrchestrator/__PrepareOutputSettings
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/orchestrators/sequential_orchestrator/SequentialOrchestrator/__GetStagesCheckpointFolder
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def Run(self) -> None:`

# Calls

- [__GetStagesToCheckpointList](../../../../../functions/python_scripts/orchestrators/sequential_orchestrator/SequentialOrchestrator/__GetStagesToCheckpointList.md)
- [GetProject](../../../../../functions/python_scripts/orchestrators/orchestrator/Orchestrator/GetProject.md)
- [IsString](../../../../../functions/sources/kratos_parameters/Parameters/IsString.md)
- [GetString](../../../../../functions/sources/kratos_parameters/Parameters/GetString.md)
- [__GetExecutionList](../../../../../functions/python_scripts/orchestrators/sequential_orchestrator/SequentialOrchestrator/__GetExecutionList.md)
- [CheckStageSettings](../../../../../functions/python_scripts/orchestrators/orchestrator/Orchestrator/CheckStageSettings.md)
- [CreateStage](../../../../../functions/python_scripts/orchestrators/orchestrator/Orchestrator/CreateStage.md)
- [RunCurrentStagePreprocess](../../../../../functions/python_scripts/orchestrators/orchestrator/Orchestrator/RunCurrentStagePreprocess.md)
- [GetOutputData](../../../../../functions/python_scripts/project/Project/GetOutputData.md)
- [GetFinalData](../../../../../functions/python_scripts/analysis_stage/AnalysisStage/GetFinalData.md)
- [RunCurrentStagePostprocess](../../../../../functions/python_scripts/orchestrators/orchestrator/Orchestrator/RunCurrentStagePostprocess.md)
- [__PrepareOutputSettings](../../../../../functions/python_scripts/orchestrators/sequential_orchestrator/SequentialOrchestrator/__PrepareOutputSettings.md)
- [__GetStagesCheckpointFolder](../../../../../functions/python_scripts/orchestrators/sequential_orchestrator/SequentialOrchestrator/__GetStagesCheckpointFolder.md)