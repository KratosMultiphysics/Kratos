---
type: C++ Method
title: RegisterVariables
resource: sources/variables.cpp#L478-L952
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/deprecated_variables/KratosApplication/RegisterDeprecatedVariables
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/cfd_variables/KratosApplication/RegisterCFDVariables
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/mesh_moving_variables/KratosApplication/RegisterALEVariables
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/mapping_variables/KratosApplication/RegisterMappingVariables
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/dem_variables/KratosApplication/RegisterDEMVariables
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/fsi_variables/KratosApplication/RegisterFSIVariables
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/mat_variables/KratosApplication/RegisterMATVariables
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/sources/kratos_application/KratosApplication/RegisterKratosCore
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void KratosApplication::RegisterVariables()`

# Calls

- [RegisterDeprecatedVariables](../../../../functions/sources/deprecated_variables/KratosApplication/RegisterDeprecatedVariables.md)
- [RegisterCFDVariables](../../../../functions/sources/cfd_variables/KratosApplication/RegisterCFDVariables.md)
- [RegisterALEVariables](../../../../functions/sources/mesh_moving_variables/KratosApplication/RegisterALEVariables.md)
- [RegisterMappingVariables](../../../../functions/sources/mapping_variables/KratosApplication/RegisterMappingVariables.md)
- [RegisterDEMVariables](../../../../functions/sources/dem_variables/KratosApplication/RegisterDEMVariables.md)
- [RegisterFSIVariables](../../../../functions/sources/fsi_variables/KratosApplication/RegisterFSIVariables.md)
- [RegisterMATVariables](../../../../functions/sources/mat_variables/KratosApplication/RegisterMATVariables.md)

# Called by

- [RegisterKratosCore](../../../../functions/sources/kratos_application/KratosApplication/RegisterKratosCore.md)