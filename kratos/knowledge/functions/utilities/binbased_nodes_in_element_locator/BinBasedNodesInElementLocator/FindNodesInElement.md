---
type: C++ Method
title: FindNodesInElement
resource: utilities/binbased_nodes_in_element_locator.h#L133-L175
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`unsigned int FindNodesInElement(Element::Pointer& pelement, DenseVector<int>& positions, Matrix& Nmat, const unsigned int max_results, PointIterator work_results, DistanceIterator work_distances, Node& work_point // is it needed from outside??????????????? )`

# Calls

- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)