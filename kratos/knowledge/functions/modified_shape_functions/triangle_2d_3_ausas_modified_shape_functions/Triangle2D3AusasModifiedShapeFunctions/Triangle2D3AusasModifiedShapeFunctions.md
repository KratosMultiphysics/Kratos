---
type: C++ Method
title: Triangle2D3AusasModifiedShapeFunctions
resource: modified_shape_functions/triangle_2d_3_ausas_modified_shape_functions.cpp#L25-L32
visibility: private
generated:
  by: okf-rs/0.7.0
---

# Signature

`Triangle2D3AusasModifiedShapeFunctions::Triangle2D3AusasModifiedShapeFunctions(const GeometryPointerType pInputGeometry, const Vector& rNodalDistances) : AusasModifiedShapeFunctions(pInputGeometry, rNodalDistances), mpTriangleSplitter(Kratos::make_shared<DivideTriangle2D3<Node>>(*pInputGeometry, rNodalDistances))`