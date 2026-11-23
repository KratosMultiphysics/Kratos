---
type: C++ Function
title: PreparePointsSearch
resource: utilities/search_utilities.h#L495-L536
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/node/SolutionStepsDataHas
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/utilities/search_utilities/PrepareSearch
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static std::vector<typename PointObject<typename TContainer::value_type>::Pointer> PreparePointsSearch( const TContainer& rStructure, const int Rank = -1 )`

# Calls

- [SolutionStepsDataHas](../../../functions/includes/node/SolutionStepsDataHas.md)

# Called by

- [PrepareSearch](../../../functions/utilities/search_utilities/PrepareSearch.md)