---
type: C++ Method
title: FindIntersectedGeometricalObjectsProcess
resource: processes/find_intersected_geometrical_objects_process.cpp#L50-L77
visibility: private
generated:
  by: okf-rs/0.7.0
---

# Signature

`FindIntersectedGeometricalObjectsProcess::FindIntersectedGeometricalObjectsProcess( Model& rModel, Parameters ThisParameters ) : mrModelPartIntersected(rModel.GetModelPart(ThisParameters["intersected_model_part_name"].GetString())), mrModelPartIntersecting(rModel.GetModelPart(ThisParameters["intersecting_model_part_name"].GetString())), mpOctree(new OctreeType())`