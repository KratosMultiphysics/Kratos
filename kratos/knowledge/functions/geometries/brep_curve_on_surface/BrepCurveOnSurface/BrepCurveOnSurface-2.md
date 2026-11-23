---
type: C++ Method
title: BrepCurveOnSurface
resource: geometries/brep_curve_on_surface.h#L99-L111
generated:
  by: okf-rs/0.7.0
---

# Signature

`BrepCurveOnSurface( typename NurbsSurfaceType::Pointer pSurface, typename NurbsCurveType::Pointer pCurve, NurbsInterval CurveNurbsInterval, bool SameCurveDirection = true) : BaseType(PointsArrayType(), &msGeometryData) , mpCurveOnSurface( Kratos::make_shared<NurbsCurveOnSurfaceType>( pSurface, pCurve)) , mCurveNurbsInterval(CurveNurbsInterval) , mSameCurveDirection(SameCurveDirection)`