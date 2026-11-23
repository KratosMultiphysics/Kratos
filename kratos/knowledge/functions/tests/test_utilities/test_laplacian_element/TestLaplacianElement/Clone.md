---
type: C++ Method
title: Clone
resource: tests/test_utilities/test_laplacian_element.cpp#L86-L97
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`Element::Pointer TestLaplacianElement::Clone( IndexType NewId, NodesArrayType const& rThisNodes ) const`

# Calls

- [GetGeometry](../../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)