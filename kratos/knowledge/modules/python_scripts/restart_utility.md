---
type: Python Module
title: restart_utility
resource: python_scripts/restart_utility.py#L1-L278
generated:
  by: okf-rs/0.7.0
relationships:
  imports:
  - target: external/kratosmultiphysics
    resolved_by: tree-sitter
    confidence: exact
  - target: external/os
    resolved_by: tree-sitter
    confidence: exact
  - target: external/pathlib
    resolved_by: tree-sitter
    confidence: exact
  member_of:
  - target: packages/KratosMultiphysics
    resolved_by: tree-sitter
    confidence: exact
---

# Contains

- [RestartUtility](../../classes/python_scripts/restart_utility/RestartUtility.md)
- [__init__](../../functions/python_scripts/restart_utility/RestartUtility/__init__.md)
- [LoadRestart](../../functions/python_scripts/restart_utility/RestartUtility/LoadRestart.md)
- [SaveRestart](../../functions/python_scripts/restart_utility/RestartUtility/SaveRestart.md)
- [IsRestartOutputStep](../../functions/python_scripts/restart_utility/RestartUtility/IsRestartOutputStep.md)
- [CreateOutputFolder](../../functions/python_scripts/restart_utility/RestartUtility/CreateOutputFolder.md)
- [_GetRestartFiles](../../functions/python_scripts/restart_utility/RestartUtility/_GetRestartFiles.md)
- [_UpdateRestartFilesMap](../../functions/python_scripts/restart_utility/RestartUtility/_UpdateRestartFilesMap.md)
- [_ClearObsoleteRestartFiles](../../functions/python_scripts/restart_utility/RestartUtility/_ClearObsoleteRestartFiles.md)
- [_GetFileLabelLoad](../../functions/python_scripts/restart_utility/RestartUtility/_GetFileLabelLoad.md)
- [_GetFileLabelSave](../../functions/python_scripts/restart_utility/RestartUtility/_GetFileLabelSave.md)
- [_ExecuteAfterLoad](../../functions/python_scripts/restart_utility/RestartUtility/_ExecuteAfterLoad.md)
- [_GetFileNamePattern](../../functions/python_scripts/restart_utility/RestartUtility/_GetFileNamePattern.md)
- [_GetSerializerFlags](../../functions/python_scripts/restart_utility/RestartUtility/_GetSerializerFlags.md)
- [_GetDefaultParameters](../../functions/python_scripts/restart_utility/RestartUtility/_GetDefaultParameters.md)
- [__GetFolderPathLoad](../../functions/python_scripts/restart_utility/RestartUtility/__GetFolderPathLoad.md)
- [__GetFolderPathSave](../../functions/python_scripts/restart_utility/RestartUtility/__GetFolderPathSave.md)
- [__GetFileNameLoad](../../functions/python_scripts/restart_utility/RestartUtility/__GetFileNameLoad.md)
- [__GetFileNameSave](../../functions/python_scripts/restart_utility/RestartUtility/__GetFileNameSave.md)
- [__GetPrettyTime](../../functions/python_scripts/restart_utility/RestartUtility/__GetPrettyTime.md)

# Imports

- `KratosMultiphysics`
- `os`
- `pathlib`

# Member of

- [KratosMultiphysics](../../packages/KratosMultiphysics.md)