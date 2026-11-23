---
type: C++ Method
title: FindIntersectedGeometricalObjectsWithOBBProcess
resource: processes/find_intersected_geometrical_objects_with_obb_process.cpp#L58-L113
visibility: private
generated:
  by: okf-rs/0.7.0
---

# Signature

`FindIntersectedGeometricalObjectsWithOBBProcess::FindIntersectedGeometricalObjectsWithOBBProcess( Model& rModel, Parameters ThisParameters ) : BaseType(rModel.GetModelPart(ThisParameters["intersected_model_part_name"].GetString()), rModel.GetModelPart(ThisParameters["intersecting_model_part_name"].GetString()))`