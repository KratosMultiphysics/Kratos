---
type: C++ Method
title: SolveIncludes
resource: sources/kratos_parameters.cpp#L244-L302
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/containers/pointer_vector_map_iterator/PointerVectorMapIterator/key
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_filesystem/FilesystemExtensions/ResolveSymlinks
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/ReadFile
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/pointer_vector_set/PointerVectorSet/pop_back
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/sources/kratos_parameters/Parameters/Parameters
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void Parameters::SolveIncludes(nlohmann::json& rJson, const std::filesystem::path& rFileName, std::vector<std::filesystem::path>& rIncludeSequence)`

# Calls

- [key](../../../../functions/containers/pointer_vector_map_iterator/PointerVectorMapIterator/key.md)
- [ResolveSymlinks](../../../../functions/sources/kratos_filesystem/FilesystemExtensions/ResolveSymlinks.md)
- [ReadFile](../../../../functions/sources/kratos_parameters/Parameters/ReadFile.md)
- [pop_back](../../../../functions/containers/pointer_vector_set/PointerVectorSet/pop_back.md)

# Called by

- [Parameters](../../../../functions/sources/kratos_parameters/Parameters/Parameters.md)