---
type: C++ Method
title: NurbsSurfaceGeometry
resource: geometries/nurbs_surface_geometry.h#L87-L109
generated:
  by: okf-rs/0.7.0
---

# Signature

`NurbsSurfaceGeometry( const PointsArrayType& rThisPoints, const SizeType PolynomialDegreeU, const SizeType PolynomialDegreeV, const Vector& rKnotsU, const Vector& rKnotsV, const Vector& rWeights, const ProjectionAlgorithm ProjectionAlgorithm = ProjectionAlgorithm::NewtonRaphson) : BaseType(rThisPoints, &msGeometryData) , mPolynomialDegreeU(PolynomialDegreeU) , mPolynomialDegreeV(PolynomialDegreeV) , mKnotsU(rKnotsU) , mKnotsV(rKnotsV) , mWeights(rWeights) , mProjectionAlgorithm(ProjectionAlgorithm)`