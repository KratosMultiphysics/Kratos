---
type: C++ Module
title: runtime_dependency_handler
resource: tests/test_utilities/runtime_dependency_handler.cpp#L1-L92
generated:
  by: okf-rs/0.7.0
relationships:
  imports:
  - target: external/runtime-dependency-handler-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-define-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/windows-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/dlfcn-h
    resolved_by: tree-sitter
    confidence: exact
  member_of:
  - target: packages/KratosMultiphysics
    resolved_by: tree-sitter
    confidence: exact
---

# Contains

- [RuntimeDependencyHandler](../../../functions/tests/test_utilities/runtime_dependency_handler/RuntimeDependencyHandler/RuntimeDependencyHandler.md)
- [LoadDependency](../../../functions/tests/test_utilities/runtime_dependency_handler/RuntimeDependencyHandler/LoadDependency.md)
- [CreateApplications](../../../functions/tests/test_utilities/runtime_dependency_handler/RuntimeDependencyHandler/CreateApplications.md)
- [ReleaseDependencies](../../../functions/tests/test_utilities/runtime_dependency_handler/RuntimeDependencyHandler/ReleaseDependencies.md)

# Imports

- `"runtime_dependency_handler.h"`
- `"includes/define.h"`
- `<windows.h>`
- `<dlfcn.h>`

# Member of

- [KratosMultiphysics](../../../packages/KratosMultiphysics.md)