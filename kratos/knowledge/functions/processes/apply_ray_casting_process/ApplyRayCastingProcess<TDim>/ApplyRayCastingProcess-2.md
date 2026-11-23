---
type: C++ Method
title: ApplyRayCastingProcess
resource: processes/apply_ray_casting_process.cpp#L43-L56
visibility: private
generated:
  by: okf-rs/0.7.0
---

# Signature

`ApplyRayCastingProcess<TDim>::ApplyRayCastingProcess( ModelPart& rVolumePart, ModelPart& rSkinPart, const double RelativeTolerance) : mRelativeTolerance(RelativeTolerance), mpFindIntersectedObjectsProcess(new FindIntersectedGeometricalObjectsProcess(rVolumePart, rSkinPart)), mIsSearchStructureAllocated(true)`