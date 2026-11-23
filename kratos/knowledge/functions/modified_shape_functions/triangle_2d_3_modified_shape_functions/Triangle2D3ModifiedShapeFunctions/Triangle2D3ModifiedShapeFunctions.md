---
type: C++ Method
title: Triangle2D3ModifiedShapeFunctions
resource: modified_shape_functions/triangle_2d_3_modified_shape_functions.cpp#L25-L32
visibility: private
generated:
  by: okf-rs/0.7.0
---

# Signature

`Triangle2D3ModifiedShapeFunctions::Triangle2D3ModifiedShapeFunctions(const GeometryPointerType pInputGeometry, const Vector& rNodalDistances) : ModifiedShapeFunctions(pInputGeometry, rNodalDistances), mpTriangleSplitter(Kratos::make_shared<DivideTriangle2D3<Node>>(*pInputGeometry, rNodalDistances))`