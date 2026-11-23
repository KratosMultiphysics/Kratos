---
type: C++ Method
title: Tetrahedra3D4AusasModifiedShapeFunctions
resource: modified_shape_functions/tetrahedra_3d_4_ausas_modified_shape_functions.cpp#L25-L32
visibility: private
generated:
  by: okf-rs/0.7.0
---

# Signature

`Tetrahedra3D4AusasModifiedShapeFunctions::Tetrahedra3D4AusasModifiedShapeFunctions(const GeometryPointerType pInputGeometry, const Vector& rNodalDistances) : AusasModifiedShapeFunctions(pInputGeometry, rNodalDistances), mpTetrahedraSplitter(Kratos::make_shared<DivideTetrahedra3D4<Node>>(*pInputGeometry, rNodalDistances))`