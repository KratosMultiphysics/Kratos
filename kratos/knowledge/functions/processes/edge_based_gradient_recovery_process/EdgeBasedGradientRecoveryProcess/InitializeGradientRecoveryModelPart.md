---
type: C++ Method
title: InitializeGradientRecoveryModelPart
resource: processes/edge_based_gradient_recovery_process.h#L451-L511
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/containers/model/Model/HasModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/model/Model/DeleteModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetBool
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/variable_utils/SetFlag
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part/ModelPart/CreateNewProperties
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/processes/edge_based_gradient_recovery_process/EdgeBasedGradientRecoveryProcess/ExecuteInitialize
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/edge_based_gradient_recovery_process/EdgeBasedGradientRecoveryProcess/InitializeGradientRecoveryStrategy
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`virtual void InitializeGradientRecoveryModelPart(ModelPart& rOriginModelPart)`

# Calls

- [HasModelPart](../../../../functions/containers/model/Model/HasModelPart.md)
- [DeleteModelPart](../../../../functions/containers/model/Model/DeleteModelPart.md)
- [GetBool](../../../../functions/sources/kratos_parameters/Parameters/GetBool.md)
- [SetFlag](../../../../functions/utilities/variable_utils/SetFlag.md)
- [CreateNewProperties](../../../../functions/sources/model_part/ModelPart/CreateNewProperties.md)

# Called by

- [ExecuteInitialize](../../../../functions/processes/edge_based_gradient_recovery_process/EdgeBasedGradientRecoveryProcess/ExecuteInitialize.md)
- [InitializeGradientRecoveryStrategy](../../../../functions/processes/edge_based_gradient_recovery_process/EdgeBasedGradientRecoveryProcess/InitializeGradientRecoveryStrategy.md)