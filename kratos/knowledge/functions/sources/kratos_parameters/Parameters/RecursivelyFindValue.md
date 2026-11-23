---
type: C++ Method
title: RecursivelyFindValue
resource: sources/kratos_parameters.cpp#L1219-L1235
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/PrettyPrintJsonString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/pointer_vector_map_iterator/PointerVectorMapIterator/key
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void Parameters::RecursivelyFindValue( const nlohmann::json& rBaseValue, const nlohmann::json& rValueToFind ) const`

# Calls

- [PrettyPrintJsonString](../../../../functions/sources/kratos_parameters/Parameters/PrettyPrintJsonString.md)
- [key](../../../../functions/containers/pointer_vector_map_iterator/PointerVectorMapIterator/key.md)