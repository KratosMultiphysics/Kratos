---
type: C++ Function
title: ReplaceEntities
resource: processes/replace_elements_and_condition_process.cpp#L36-L60
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/processes/replace_elements_and_condition_process/ReplaceElementsAndConditionsProcess/Execute
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ReplaceEntities( typename EntitiesUtilities::EntitityIdentifier<typename TContainer::value_type>& rEntitityIdentifier, TContainer& rEntityContainer, const std::unordered_set<std::size_t>& rSetOfIds )`

# Called by

- [Execute](../../../functions/processes/replace_elements_and_condition_process/ReplaceElementsAndConditionsProcess/Execute.md)