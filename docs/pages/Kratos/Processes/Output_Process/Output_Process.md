---
title: Output Process
keywords: process core Output
tags: [Output process]
sidebar: kratos_core_processes
summary: 
---

# Output Processes

This section intends to describe what is an output process in Kratos and the different interfaces that are available from python.

## Description

Output Processes in kratos, like process, are a set of utilities that will be executed during an `AnalysisStage`. The main characteristic of the output process procsses is that in addition to the regular execution points of a process, they have an additional fixed execute point `PrintOutput` and a function to evaluate if the current step is valid for printing `IsOutputStep` . For more information see [AnalysisStage Sequence Diagram](../Sequence_Diagrams/General/AnalysisStage):

- `ExecuteInitialize`: Will be called during the initialize sequence of an `AnalysisStage`, before the initialization of the `Solver`.

- `ExecuteBeforeSolustionLoop`: Will be called during the initialize sequence of an `AnalysisStage`, after the initialization of the `Solver`

- `ExecuteInitializeSolutionStep`: Will be called at the begining of each solution loop, before executing the preconditioners and solvers

- `ExecuteFinalizeSolutionStep`: Will be called at the end of each solution loop, after executing the preconditioners and solvers but before the output stage.

- `ExecuteBeforeOutputStep`: Will be called at the begining of the output stage for every output process active, before printing the results

- `ExecuteAfterOutputStep`: Will be called at the end of the output stage for every output process active, after printing the results

- `ExecuteFinalize`: Will be called during the finalize sequence of an `AnalysisStage`, just before existing the stage.

- `PrintOutput` Will ve called during the print sequence of the `AnalysisStage` if the print condition is meet.

- `IsOutputStep` Will enable (`true`) or disable (`false`) the print condition.

## List of Output Processes:

- [GiD Output Process](/GiD_Output_process)
- [Vtk Output Process](/VTK_Output_process)