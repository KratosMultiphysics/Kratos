---
type: C++ Method
title: Tetrahedra3D4ModifiedShapeFunctions
resource: modified_shape_functions/tetrahedra_3d_4_modified_shape_functions.cpp#L25-L32
visibility: private
generated:
  by: okf-rs/0.7.0
---

# Signature

`Tetrahedra3D4ModifiedShapeFunctions::Tetrahedra3D4ModifiedShapeFunctions(const GeometryPointerType pInputGeometry, const Vector& rNodalDistances) : ModifiedShapeFunctions(pInputGeometry, rNodalDistances), mpTetrahedraSplitter(Kratos::make_shared<DivideTetrahedra3D4<Node>>(*pInputGeometry, rNodalDistances))`