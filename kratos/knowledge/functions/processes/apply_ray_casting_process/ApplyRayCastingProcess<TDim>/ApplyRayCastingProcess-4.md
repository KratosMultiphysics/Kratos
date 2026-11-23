---
type: C++ Method
title: ApplyRayCastingProcess
resource: processes/apply_ray_casting_process.cpp#L73-L98
visibility: private
generated:
  by: okf-rs/0.7.0
---

# Signature

`ApplyRayCastingProcess<TDim>::ApplyRayCastingProcess( FindIntersectedGeometricalObjectsProcess &TheFindIntersectedObjectsProcess, const double RelativeTolerance, const Variable<double>* pDistanceVariable, const DistanceDatabase& rDistanceDatabase) : mRelativeTolerance(RelativeTolerance), mpFindIntersectedObjectsProcess(&TheFindIntersectedObjectsProcess), mIsSearchStructureAllocated(false), mpDistanceVariable(pDistanceVariable)`