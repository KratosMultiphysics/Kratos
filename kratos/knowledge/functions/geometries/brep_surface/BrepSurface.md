---
type: C++ Function
title: BrepSurface
resource: geometries/brep_surface.h#L145-L154
generated:
  by: okf-rs/0.7.0
---

# Signature

`explicit BrepSurface( BrepSurface<TOtherContainerPointType, TShiftedBoundary, TOtherContainerPointEmbeddedType> const& rOther ) : BaseType( rOther ) , mpNurbsSurface(rOther.mpNurbsSurface) , mOuterLoopArray(rOther.mOuterLoopArray) , mInnerLoopArray(rOther.mInnerLoopArray) , mEmbeddedEdgesArray(rOther.mEmbeddedEdgesArray) , mIsTrimmed(rOther.mIsTrimmed)`