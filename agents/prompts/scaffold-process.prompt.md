---
description: "Scaffold a new Kratos Process class for an application — generates C++ header, source, and pybind11 binding stub."
name: "Scaffold Kratos Process"
argument-hint: "Process name and target application, e.g. 'ComputeVolumetricStrainProcess in StructuralMechanicsApplication'"
agent: agent
tools: [read, edit, search]
---

Generate a complete scaffold for a new Kratos `Process` in this repository.

## Inputs

The user has provided: **{{input}}**

Parse the following from that input:
- **ProcessName** — PascalCase name of the process (e.g. `ComputeVolumetricStrainProcess`)
- **ApplicationName** — the target application folder under `applications/` (e.g. `StructuralMechanicsApplication`)

## Steps

1. **Read conventions** — load [cpp-conventions instructions](../instructions/cpp-conventions.instructions.md) and [python-conventions instructions](../instructions/python-conventions.instructions.md).

2. **Explore the target application** — search for an existing process in `applications/<ApplicationName>/custom_processes/` to use as a style reference.

3. **Generate the header** — create `applications/<ApplicationName>/custom_processes/<process_name>.h` with:
   - `#pragma once`
   - Correct namespace (`Kratos`)
   - Class inheriting from `Process`
   - `KRATOS_CLASS_POINTER_DEFINITION` macro
   - Constructor accepting `ModelPart&` and `Parameters`
   - Declarations for `ExecuteInitialize()`, `ExecuteInitializeSolutionStep()`, `Execute()`, `ExecuteFinalizeSolutionStep()`, `ExecuteFinalize()`

4. **Generate the source** — create `applications/<ApplicationName>/custom_processes/<process_name>.cpp` with:
   - `KRATOS_TRY` / `KRATOS_CATCH("")` in every method body
   - `KRATOS_INFO` log at start of `Execute()`
   - Member storage of `mrModelPart` and `mParameters`

5. **Generate the binding stub** — append the new process to `applications/<ApplicationName>/custom_python/add_custom_processes_to_python.cpp` following existing entries:
   - `py::class_<ProcessName, ProcessName::Pointer, Process>(m, "ProcessName")`
   - Expose constructor and all stage methods with docstrings

6. **Summarise** what was created and remind the user to:
   - Register the process in `<application_name>_application.cpp` if required
   - Add the `.cpp` source to `CMakeLists.txt` if not covered by glob
