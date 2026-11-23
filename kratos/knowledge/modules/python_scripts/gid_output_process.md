---
type: Python Module
title: gid_output_process
resource: python_scripts/gid_output_process.py#L1-L746
generated:
  by: okf-rs/0.7.0
relationships:
  imports:
  - target: external/kratosmultiphysics-as-km
    resolved_by: tree-sitter
    confidence: exact
  - target: external/kratosmultiphysics
    resolved_by: tree-sitter
    confidence: exact
  - target: external/kratosmultiphysics-deprecation-management
    resolved_by: tree-sitter
    confidence: exact
  - target: external/os
    resolved_by: tree-sitter
    confidence: exact
  - target: external/pathlib
    resolved_by: tree-sitter
    confidence: exact
  - target: external/kratosmultiphysics-mpi-distributed-gid-output-process
    resolved_by: tree-sitter
    confidence: exact
  - target: external/kratosmultiphysics-point-output-process
    resolved_by: tree-sitter
    confidence: exact
  member_of:
  - target: packages/KratosMultiphysics
    resolved_by: tree-sitter
    confidence: exact
---

# Contains

- [Factory](../../functions/python_scripts/gid_output_process/Factory.md)
- [GiDOutputProcess](../../classes/python_scripts/gid_output_process/GiDOutputProcess.md)
- [__init__](../../functions/python_scripts/gid_output_process/GiDOutputProcess/__init__.md)
- [TranslateLegacyVariablesAccordingToCurrentStandard](../../functions/python_scripts/gid_output_process/GiDOutputProcess/TranslateLegacyVariablesAccordingToCurrentStandard.md)
- [ExecuteInitialize](../../functions/python_scripts/gid_output_process/GiDOutputProcess/ExecuteInitialize.md)
- [ExecuteBeforeSolutionLoop](../../functions/python_scripts/gid_output_process/GiDOutputProcess/ExecuteBeforeSolutionLoop.md)
- [ExecuteInitializeSolutionStep](../../functions/python_scripts/gid_output_process/GiDOutputProcess/ExecuteInitializeSolutionStep.md)
- [ExecuteFinalizeSolutionStep](../../functions/python_scripts/gid_output_process/GiDOutputProcess/ExecuteFinalizeSolutionStep.md)
- [IsOutputStep](../../functions/python_scripts/gid_output_process/GiDOutputProcess/IsOutputStep.md)
- [PrintOutput](../../functions/python_scripts/gid_output_process/GiDOutputProcess/PrintOutput.md)
- [ExecuteFinalize](../../functions/python_scripts/gid_output_process/GiDOutputProcess/ExecuteFinalize.md)
- [_InitializeGiDIO](../../functions/python_scripts/gid_output_process/GiDOutputProcess/_InitializeGiDIO.md)
- [__get_pretty_time](../../functions/python_scripts/gid_output_process/GiDOutputProcess/__get_pretty_time.md)
- [__get_gidpost_flag](../../functions/python_scripts/gid_output_process/GiDOutputProcess/__get_gidpost_flag.md)
- [_initialize_cut_output](../../functions/python_scripts/gid_output_process/GiDOutputProcess/_initialize_cut_output.md)
- [_InitializeListFiles](../../functions/python_scripts/gid_output_process/GiDOutputProcess/_InitializeListFiles.md)
- [__define_output_plane](../../functions/python_scripts/gid_output_process/GiDOutputProcess/__define_output_plane.md)
- [__write_mesh](../../functions/python_scripts/gid_output_process/GiDOutputProcess/__write_mesh.md)
- [__initialize_results](../../functions/python_scripts/gid_output_process/GiDOutputProcess/__initialize_results.md)
- [__write_nodal_results](../../functions/python_scripts/gid_output_process/GiDOutputProcess/__write_nodal_results.md)
- [__write_gp_results](../../functions/python_scripts/gid_output_process/GiDOutputProcess/__write_gp_results.md)
- [__write_nonhistorical_nodal_results](../../functions/python_scripts/gid_output_process/GiDOutputProcess/__write_nonhistorical_nodal_results.md)
- [__write_nodal_flags](../../functions/python_scripts/gid_output_process/GiDOutputProcess/__write_nodal_flags.md)
- [__write_elemental_conditional_flags](../../functions/python_scripts/gid_output_process/GiDOutputProcess/__write_elemental_conditional_flags.md)
- [__finalize_results](../../functions/python_scripts/gid_output_process/GiDOutputProcess/__finalize_results.md)
- [__write_step_to_list](../../functions/python_scripts/gid_output_process/GiDOutputProcess/__write_step_to_list.md)
- [__restart_list_files](../../functions/python_scripts/gid_output_process/GiDOutputProcess/__restart_list_files.md)
- [__remove_list_files](../../functions/python_scripts/gid_output_process/GiDOutputProcess/__remove_list_files.md)
- [__remove_post_results_files](../../functions/python_scripts/gid_output_process/GiDOutputProcess/__remove_post_results_files.md)
- [_CreateCuttingUtility](../../functions/python_scripts/gid_output_process/GiDOutputProcess/_CreateCuttingUtility.md)
- [_SetCurrentTimeParameters](../../functions/python_scripts/gid_output_process/GiDOutputProcess/_SetCurrentTimeParameters.md)
- [_ValidateFileName](../../functions/python_scripts/gid_output_process/GiDOutputProcess/_ValidateFileName.md)

# Imports

- `KratosMultiphysics as KM`
- `KratosMultiphysics`
- `KratosMultiphysics.deprecation_management`
- `os`
- `pathlib`
- `KratosMultiphysics.mpi.distributed_gid_output_process`
- `KratosMultiphysics.point_output_process`

# Member of

- [KratosMultiphysics](../../packages/KratosMultiphysics.md)