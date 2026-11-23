---
type: C++ Method
title: QuadraturePointGeometry
resource: geometries/quadrature_point_geometry.h#L157-L173
generated:
  by: okf-rs/0.7.0
---

# Signature

`QuadraturePointGeometry( const PointsArrayType& ThisPoints, const IntegrationPointType& ThisIntegrationPoint, const Matrix& ThisShapeFunctionsValues, const DenseVector<Matrix>& ThisShapeFunctionsDerivatives, GeometryType* pGeometryParent) : BaseType(ThisPoints, &mGeometryData) , mGeometryData( &msGeometryDimension, GeometryShapeFunctionContainerType( GeometryData::IntegrationMethod::GI_GAUSS_1, ThisIntegrationPoint, ThisShapeFunctionsValues, ThisShapeFunctionsDerivatives)) , mpGeometryParent(pGeometryParent)`