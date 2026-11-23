---
type: C++ Method
title: BrepCurveOnSurface
resource: geometries/brep_curve_on_surface.h#L114-L122
generated:
  by: okf-rs/0.7.0
---

# Signature

`BrepCurveOnSurface( NurbsCurveOnSurfacePointerType pNurbsCurveOnSurface, bool SameCurveDirection = true) : BaseType(PointsArrayType(), &msGeometryData) , mpCurveOnSurface(pNurbsCurveOnSurface) , mCurveNurbsInterval(pNurbsCurveOnSurface->DomainInterval()) , mSameCurveDirection(SameCurveDirection)`