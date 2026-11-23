---
type: C++ Method
title: BrepSurface
resource: geometries/brep_surface.h#L114-L125
generated:
  by: okf-rs/0.7.0
---

# Signature

`BrepSurface( typename NurbsSurfaceType::Pointer pSurface, BrepCurveOnSurfaceLoopArrayType& BrepOuterLoopArray, BrepCurveOnSurfaceLoopArrayType& BrepInnerLoopArray, bool IsTrimmed) : BaseType(PointsArrayType(), &msGeometryData) , mpNurbsSurface(pSurface) , mOuterLoopArray(BrepOuterLoopArray) , mInnerLoopArray(BrepInnerLoopArray) , mIsTrimmed(IsTrimmed)`