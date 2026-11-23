---
type: C++ Method
title: GatherModelPartUtility
resource: mpi/utilities/gather_modelpart_utility.cpp#L29-L148
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/geometry/PointsNumber
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`GatherModelPartUtility::GatherModelPartUtility( const int GatherRank, ModelPart& rOriginModelPart, const int MeshId, ModelPart& rDestinationModelPart ) : mrModelPart(rDestinationModelPart), mGatherRank(GatherRank)`

# Calls

- [GetGeometry](../../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [PointsNumber](../../../../../functions/geometries/geometry/PointsNumber.md)