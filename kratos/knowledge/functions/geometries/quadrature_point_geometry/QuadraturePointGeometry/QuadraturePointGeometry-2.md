---
type: C++ Method
title: QuadraturePointGeometry
resource: geometries/quadrature_point_geometry.h#L98-L113
generated:
  by: okf-rs/0.7.0
---

# Signature

`QuadraturePointGeometry( const PointsArrayType& ThisPoints, const IntegrationPointsContainerType& rIntegrationPoints, const ShapeFunctionsValuesContainerType& rShapeFunctionValues, const ShapeFunctionsLocalGradientsContainerType& rShapeFunctionsDerivativesVector, GeometryType* pGeometryParent) : BaseType(ThisPoints, &mGeometryData) , mGeometryData( &msGeometryDimension, GeometryData::IntegrationMethod::GI_GAUSS_1, rIntegrationPoints, rShapeFunctionValues, rShapeFunctionsDerivativesVector) , mpGeometryParent(pGeometryParent)`