---
type: C++ Method
title: SplitBoundary
resource: processes/split_internal_interfaces_process.cpp#L71-L188
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/geometry/GenerateBoundariesEntities
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part/ModelPart/CreateNewProperties
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/geometry/PointsNumber
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/processes/split_internal_interfaces_process/SplitInternalInterfacesProcess/ExecuteInitialize
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void SplitInternalInterfacesProcess::SplitBoundary( const std::size_t PropertyIdBeingProcessed, const std::size_t InterfaceConditionsPropertyId, const std::unordered_map<IndexType, IndexType>& rConditionsParentMap, ModelPart& rModelPart)`

# Calls

- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [GenerateBoundariesEntities](../../../../functions/geometries/geometry/GenerateBoundariesEntities.md)
- [CreateNewProperties](../../../../functions/sources/model_part/ModelPart/CreateNewProperties.md)
- [PointsNumber](../../../../functions/geometries/geometry/PointsNumber.md)

# Called by

- [ExecuteInitialize](../../../../functions/processes/split_internal_interfaces_process/SplitInternalInterfacesProcess/ExecuteInitialize.md)