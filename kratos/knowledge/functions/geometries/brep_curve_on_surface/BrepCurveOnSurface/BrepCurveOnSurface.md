---
type: C++ Method
title: BrepCurveOnSurface
resource: geometries/brep_curve_on_surface.h#L85-L96
generated:
  by: okf-rs/0.7.0
---

# Signature

`BrepCurveOnSurface( typename NurbsSurfaceType::Pointer pSurface, typename NurbsCurveType::Pointer pCurve, bool SameCurveDirection = true) : BaseType(PointsArrayType(), &msGeometryData) , mpCurveOnSurface( Kratos::make_shared<NurbsCurveOnSurfaceType>( pSurface, pCurve)) , mCurveNurbsInterval(pCurve->DomainInterval()) , mSameCurveDirection(SameCurveDirection)`