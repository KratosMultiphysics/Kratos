---
type: C++ Method
title: ReadSolid
resource: input_output/stl_io.cpp#L370-L413
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/stl_io/StlIO/ReadFacet
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/input_output/stl_io/StlIO/ReadModelPart
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void StlIO::ReadSolid( ModelPart& rThisModelPart, const std::function<void(ModelPart&, NodesArrayType&)>& rCreateEntityFunctor )`

# Calls

- [ReadFacet](../../../../functions/input_output/stl_io/StlIO/ReadFacet.md)

# Called by

- [ReadModelPart](../../../../functions/input_output/stl_io/StlIO/ReadModelPart.md)