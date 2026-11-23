---
type: C++ Method
title: BrepVolume
resource: geometries/brep_volume.h#L107-L118
generated:
  by: okf-rs/0.7.0
---

# Signature

`BrepVolume( typename NurbsVolumeType::Pointer pVolume, BrepSurfaceOnVolumeLoopArrayType& BrepOuterLoopArray, BrepSurfaceOnVolumeLoopArrayType& BrepInnerLoopArray, bool IsTrimmed) : BaseType(PointsArrayType(), &msGeometryData) , mpNurbsVolume(pVolume) , mOuterLoopArray(BrepOuterLoopArray) , mInnerLoopArray(BrepInnerLoopArray) , mIsTrimmed(IsTrimmed)`