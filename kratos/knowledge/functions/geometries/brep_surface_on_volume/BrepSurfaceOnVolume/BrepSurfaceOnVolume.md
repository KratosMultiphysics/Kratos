---
type: C++ Method
title: BrepSurfaceOnVolume
resource: geometries/brep_surface_on_volume.h#L88-L98
generated:
  by: okf-rs/0.7.0
---

# Signature

`BrepSurfaceOnVolume( typename NurbsVolumeType::Pointer pVolume, typename NurbsSurfaceType::Pointer pSurface) : BaseType(PointsArrayType(), &msGeometryData) , mpSurfaceOnVolume( Kratos::make_shared<NurbsSurfaceOnVolumeType>( pVolume, pSurface)) , mSurfaceNurbsIntervalU(pSurface->DomainIntervalU()) , mSurfaceNurbsIntervalV(pSurface->DomainIntervalV())`