---
type: C++ Method
title: SetUpSystemIds
resource: future/solving_strategies/schemes/implicit_scheme.h#L328-L347
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/dof_utilities/dof_array_utilities/DofArrayUtilities/SetDofEquationIds
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/dof_utilities/dof_array_utilities/DofArrayUtilities/SetEffectiveDofEquationIds
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/future/solving_strategies/schemes/implicit_scheme/ImplicitScheme/InitializeLinearSystem
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`virtual void SetUpSystemIds( typename DofsArrayType::Pointer pDofSet, typename DofsArrayType::Pointer pEffectiveDofSet)`

# Calls

- [SetDofEquationIds](../../../../../../functions/utilities/dof_utilities/dof_array_utilities/DofArrayUtilities/SetDofEquationIds.md)
- [SetEffectiveDofEquationIds](../../../../../../functions/utilities/dof_utilities/dof_array_utilities/DofArrayUtilities/SetEffectiveDofEquationIds.md)

# Called by

- [InitializeLinearSystem](../../../../../../functions/future/solving_strategies/schemes/implicit_scheme/ImplicitScheme/InitializeLinearSystem.md)