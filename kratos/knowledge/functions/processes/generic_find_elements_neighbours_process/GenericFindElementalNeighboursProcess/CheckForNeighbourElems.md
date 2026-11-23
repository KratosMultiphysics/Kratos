---
type: C++ Method
title: CheckForNeighbourElems
resource: processes/generic_find_elements_neighbours_process.cpp#L69-L106
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/processes/generic_find_elements_neighbours_process/GenericFindElementalNeighboursProcess/Execute
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`GlobalPointer<Element> GenericFindElementalNeighboursProcess::CheckForNeighbourElems (const Geometry<Node >& rBoundaryGeom, Element & rElement, const int CurrentRank)`

# Called by

- [Execute](../../../../functions/processes/generic_find_elements_neighbours_process/GenericFindElementalNeighboursProcess/Execute.md)