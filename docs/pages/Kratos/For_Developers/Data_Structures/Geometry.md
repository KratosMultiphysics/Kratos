---
title: Geometry
keywords: 
tags: [Geometry Use]
sidebar: kratos_for_developers
summary: 
---

The integration of FE functions, typically requires the evaluation of functions at given integration points as well as the computation of the shape functions and of their derivates on such points

the following code snippet is useful to compute some elementary FE operations, abstracting the computations from the actual choice of the geometry.
The idea is that in general it is possible to write an element in a similar way for, say, a triangle and a quadrilateral.

If we are within an element, we can obtain a reference to the geometry (read as topology if you like) of the element by doing

```cpp
const GeometryType& rGeom = this->GetGeometry();
```

The geometry provides access to some very useful methods to compute the `Area`, `Volume`, etc.

It also grants access to the nodes of the element. The first node can be accessed via:

```cpp
rGeom[0]
```

or

```cpp
this->GetGeometry()[0]
```

In the same way the rest of the nodes can be accessed by using [1], [2], ...

The geometry also provides facilities to compute Jacobians, determinants, shape functions, etc. 
For example to compute the Determinant of the Jacobian and the shape function derivatives at once for all of the gauss points of interest, one can do

```cpp
Vector DetJ;
ShapeFunctionDerivativesArrayType DN_DX;
rGeom.ShapeFunctionsIntegrationPointsGradients(DN_DX,DetJ,GeometryData::GI_GAUSS_2);
```

The vector `DetJ` will contain the determinant of the Jacobian for all of the integration points considered in the integration rule considered.

The choice of the `integration rule` is controlled here by the parameter `GeometryData::GI_GAUSS_2` which here corresponds to using the 2nd order integration rule,
which typically implies 2 gauss points per direction if kroneker product can be used, or to a similar precision for the case of triangles and tetrahedra.
GI_GAUSS_1,2,3,4 are available, providing an increasing precision at the price of using more integration points

`DN_DX[IntPointNumber](I,k)` contains the derivative of the shape function for node I, with respect to the "k" component, evaluated at the gauss point `IntPointNumber`
   
One may obtain the Shape functions at the gauss points by using:

```cpp
Matrix NContainer = rGeom.ShapeFunctionsValues(GeometryData::GI_GAUSS_2);
```

Ncontainer is defined such that `Ncontainer(IntPointNumber, I)` contains the shape functions of node I evaluated at the integration point `IntPointNumber`

Note that if the parameter `GeometryData::GI_GAUSS_2` can be omitted in all of the functions. In this case the default integration rule will be used. Such integration rule will provide exact results for the computation of a laplacian or of a the stiffness matrix for a linear structural problem.

The list of integration points can be queried by:

```cpp
const GeometryType::IntegrationPointsArrayType& IntegrationPoints = 
rGeom.IntegrationPoints(GeometryData::GI_GAUSS_2);
```

This can be used for example to compute the total area of the element as

```cpp
double area = 0.0;
for (unsigned int g = 0; g < IntegrationPoints.size(); g++)
    area += DetJ[g] * IntegrationPoints[g].Weight();
```

Similarly the geometry allows the computation of the Jacobians on all of the gauss points as well as of the "Local" ShapeFunctionDerivatives, expressed in terms of parent coordinates.
