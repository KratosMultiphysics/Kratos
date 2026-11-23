---
type: C++ Function
title: AddDofWithReaction
resource: utilities/variable_utils.h#L1403-L1429
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/node/SolutionStepsDataHas
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/variable_utils/CheckVariableExists
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void AddDofWithReaction( const TVarType& rVar, const TVarType& rReactionVar, ModelPart& rModelPart )`

# Calls

- [SolutionStepsDataHas](../../../functions/includes/node/SolutionStepsDataHas.md)
- [CheckVariableExists](../../../functions/utilities/variable_utils/CheckVariableExists.md)