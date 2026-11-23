---
type: C++ Method
title: SetUpDofArrays
resource: future/solving_strategies/schemes/implicit_scheme.h#L306-L320
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/dof_utilities/dof_array_utilities/DofArrayUtilities/SetUpEffectiveDofArray
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/future/solving_strategies/schemes/implicit_scheme/ImplicitScheme/InitializeLinearSystem
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`virtual std::pair<std::size_t, std::size_t> SetUpDofArrays( typename DofsArrayType::Pointer pDofSet, typename DofsArrayType::Pointer pEffectiveDofSet)`

# Calls

- [SetUpEffectiveDofArray](../../../../../../functions/utilities/dof_utilities/dof_array_utilities/DofArrayUtilities/SetUpEffectiveDofArray.md)

# Called by

- [InitializeLinearSystem](../../../../../../functions/future/solving_strategies/schemes/implicit_scheme/ImplicitScheme/InitializeLinearSystem.md)