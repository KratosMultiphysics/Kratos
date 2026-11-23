---
type: C++ Module
title: kratos_filesystem
resource: sources/kratos_filesystem.cpp#L1-L153
generated:
  by: okf-rs/0.7.0
relationships:
  imports:
  - target: external/algorithm
    resolved_by: tree-sitter
    confidence: exact
  - target: external/chrono
    resolved_by: tree-sitter
    confidence: exact
  - target: external/set
    resolved_by: tree-sitter
    confidence: exact
  - target: external/thread
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-kratos-filesystem-h
    resolved_by: tree-sitter
    confidence: exact
  member_of:
  - target: packages/KratosMultiphysics
    resolved_by: tree-sitter
    confidence: exact
---

# Contains

- [exists](../../functions/sources/kratos_filesystem/exists.md)
- [is_regular_file](../../functions/sources/kratos_filesystem/is_regular_file.md)
- [is_directory](../../functions/sources/kratos_filesystem/is_directory.md)
- [create_directory](../../functions/sources/kratos_filesystem/create_directory.md)
- [create_directories](../../functions/sources/kratos_filesystem/create_directories.md)
- [remove](../../functions/sources/kratos_filesystem/remove.md)
- [remove_all](../../functions/sources/kratos_filesystem/remove_all.md)
- [rename](../../functions/sources/kratos_filesystem/rename.md)
- [parent_path](../../functions/sources/kratos_filesystem/parent_path.md)
- [filename](../../functions/sources/kratos_filesystem/filename.md)
- [CurrentWorkingDirectory](../../functions/sources/kratos_filesystem/FilesystemExtensions/CurrentWorkingDirectory.md)
- [JoinPaths](../../functions/sources/kratos_filesystem/FilesystemExtensions/JoinPaths.md)
- [ListDirectory](../../functions/sources/kratos_filesystem/FilesystemExtensions/ListDirectory.md)
- [MPISafeCreateDirectories](../../functions/sources/kratos_filesystem/FilesystemExtensions/MPISafeCreateDirectories.md)
- [ResolveSymlinks](../../functions/sources/kratos_filesystem/FilesystemExtensions/ResolveSymlinks.md)

# Imports

- `<algorithm>`
- `<chrono>`
- `<set>`
- `<thread>`
- `"includes/kratos_filesystem.h"`

# Member of

- [KratosMultiphysics](../../packages/KratosMultiphysics.md)