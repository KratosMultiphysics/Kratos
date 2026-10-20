---
title: Checkpoint Process
keywords: checkpoint process core
tags: [checkpoint process]
sidebar: kratos_core_processes
summary: Save or load a model part checkpoint using a file serializer.
---

# Checkpoint Process

## Description

This process saves or loads the selected model part using a `FileSerializer`.
Set `checkpoint_state` to `"save"` to write a checkpoint, or to `"load"` to
restore one.

## Execution

- With `"save"`, the checkpoint is written in `ExecuteFinalizeSolutionStep`.
- With `"load"`, the checkpoint is restored in `ExecuteInitializeSolutionStep`.

The serializer key is `"check_point"`. In `file_name`, the `<time>` and
`<step>` placeholders are replaced with the model part's current `TIME` and
`STEP` values from `ProcessInfo`. Use the same resulting filename and serializer
key when saving and loading a checkpoint.

## Parameters & Defaults

```json
{
    "model_part_name" : "please_specify_model_part_name",
    "checkpoint_state": "",
    "file_name"       : ""
}
```

- **`model_part_name`**: Name of the model part whose state is saved or loaded.
- **`checkpoint_state`**: Required operation: `"save"` or `"load"`.
- **`file_name`**: Checkpoint file name. It may include `<time>` and/or
  `<step>` placeholders to create distinct names for different solution steps.

## Example

The following process settings save the model part at the end of each solution
step. To restore a checkpoint, set `checkpoint_state` to `"load"` and configure
the matching file name; loading occurs at the beginning of the solution step.

```json
{
    "python_module" : "checkpoint_process",
    "kratos_module" : "KratosMultiphysics",
    "process_name"  : "CheckPointProcess",
    "Parameters"    : {
        "model_part_name" : "MainModelPart",
        "checkpoint_state": "save",
        "file_name"       : "checkpoint_<time>_<step>"
    }
}
```
