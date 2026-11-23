---
type: C++ Method
title: CreateTetrahedraInCell
resource: modeler/entity_generation/generate_tetrahedral_elements_with_cell_color.cpp#L68-L96
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/modeler/entity_generation/generate_tetrahedral_elements_with_cell_color/GenerateTetrahedralElementsWithCellColor/Generate
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void GenerateTetrahedralElementsWithCellColor::CreateTetrahedraInCell( ModelPart::ElementsContainerType& rElements, Element::NodesArrayType& rCellNodes, const std::size_t StartId, Properties::Pointer& pProperties, const Element& rPrototype)`

# Called by

- [Generate](../../../../../functions/modeler/entity_generation/generate_tetrahedral_elements_with_cell_color/GenerateTetrahedralElementsWithCellColor/Generate.md)