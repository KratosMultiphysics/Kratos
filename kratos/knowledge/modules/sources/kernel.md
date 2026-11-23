---
type: C++ Module
title: kernel
resource: sources/kernel.cpp#L1-L269
generated:
  by: okf-rs/0.7.0
relationships:
  imports:
  - target: external/iostream
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-kernel-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-kratos-version-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-data-communicator-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-parallel-environment-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-registry-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/input-output-logger-h
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

- [Kernel](../../functions/sources/kernel/Kernel/Kernel.md)
- [Kernel](../../functions/sources/kernel/Kernel/Kernel-2.md)
- [PrintInfo](../../functions/sources/kernel/Kernel/PrintInfo.md)
- [Initialize](../../functions/sources/kernel/Kernel/Initialize.md)
- [GetLibraryList](../../functions/sources/kernel/Kernel/GetLibraryList.md)
- [IsImported](../../functions/sources/kernel/Kernel/IsImported.md)
- [IsLibraryAvailable](../../functions/sources/kernel/Kernel/IsLibraryAvailable.md)
- [IsDistributedRun](../../functions/sources/kernel/Kernel/IsDistributedRun.md)
- [ImportApplication](../../functions/sources/kernel/Kernel/ImportApplication.md)
- [Info](../../functions/sources/kernel/Kernel/Info.md)
- [PrintInfo](../../functions/sources/kernel/Kernel/PrintInfo-2.md)
- [PrintData](../../functions/sources/kernel/Kernel/PrintData.md)
- [BuildType](../../functions/sources/kernel/Kernel/BuildType.md)
- [Version](../../functions/sources/kernel/Kernel/Version.md)
- [OSName](../../functions/sources/kernel/Kernel/OSName.md)
- [PythonVersion](../../functions/sources/kernel/Kernel/PythonVersion.md)
- [Compiler](../../functions/sources/kernel/Kernel/Compiler.md)
- [SetPythonVersion](../../functions/sources/kernel/Kernel/SetPythonVersion.md)
- [PrintParallelismSupportInfo](../../functions/sources/kernel/Kernel/PrintParallelismSupportInfo.md)

# Imports

- `<iostream>`
- `"includes/kernel.h"`
- `"includes/kratos_version.h"`
- `"includes/data_communicator.h"`
- `"includes/parallel_environment.h"`
- `"includes/registry.h"`
- `"input_output/logger.h"`
- `"utilities/parallel_utilities.h"`

# Member of

- [KratosMultiphysics](../../packages/KratosMultiphysics.md)