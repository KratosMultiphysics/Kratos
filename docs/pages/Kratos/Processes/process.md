---
title: 
keywords: process core
tags: [process]
sidebar: kratos_core_processes
summary: 
---

# Processes

This section intends to describe what is a process in Kratos and the different interfaces that are available from python.

## Description

Processes in kratos are a set of utilities that will be executed during an `AnalysisStage`. The main characteristic of the procsses is that they have a set of fixed execute points (see [AnalysisStage Sequence Diagram](../Sequence_Diagrams/General/AnalysisStage)):

- `ExecuteInitialize`: Will be called during the initialize sequence of an `AnalysisStage`, before the initialization of the `Solver`.

- `ExecuteBeforeSolustionLoop`: Will be called during the initialize sequence of an `AnalysisStage`, after the initialization of the `Solver`

- `ExecuteInitializeSolutionStep`: Will be called at the begining of each solution loop, before executing the preconditioners and solvers

- `ExecuteFinalizeSolutionStep`: Will be called at the end of each solution loop, after executing the preconditioners and solvers but before the output stage.

- `ExecuteBeforeOutputStep`: Will be called at the begining of the output stage for every output process active, before printing the results

- `ExecuteAfterOutputStep`: Will be called at the end of the output stage for every output process active, after printing the results

- `ExecuteFinalize`: Will be called during the finalize sequence of an `AnalysisStage`, just before existing the stage.

## List of processes:

### Value Assignement

- [Assign Scalar Variable](/Assign_Values/assign_scalar_variable_process)
- [Assign Vector Variable](/Assign_Values/assign_vector_variable_process)
- [Assign Vector By Direction](/Assign_Values/assign_vector_by_direction_process)
- [Assign Scalar Input](/Assign_Values/assign_scalar_input_process)

### Calculation Utilities