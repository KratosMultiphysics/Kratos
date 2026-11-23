---
type: C++ Module
title: add_geometries_to_python
resource: python/add_geometries_to_python.cpp#L1-L510
generated:
  by: okf-rs/0.7.0
relationships:
  imports:
  - target: external/includes-define-python-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/geometries-point-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-node-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/python-containers-interface-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/python-add-geometries-to-python-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/geometries-geometry-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/geometries-line-2d-2-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/geometries-line-2d-3-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/geometries-triangle-2d-3-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/geometries-triangle-2d-6-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/geometries-quadrilateral-2d-4-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/geometries-quadrilateral-2d-8-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/geometries-quadrilateral-2d-9-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/geometries-line-3d-2-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/geometries-line-3d-3-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/geometries-triangle-3d-3-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/geometries-triangle-3d-6-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/geometries-quadrilateral-3d-4-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/geometries-quadrilateral-3d-8-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/geometries-quadrilateral-3d-9-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/geometries-tetrahedra-3d-4-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/geometries-tetrahedra-3d-10-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/geometries-pyramid-3d-5-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/geometries-pyramid-3d-13-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/geometries-prism-3d-6-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/geometries-prism-3d-15-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/geometries-hexahedra-3d-8-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/geometries-hexahedra-3d-20-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/geometries-hexahedra-3d-27-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/geometries-nurbs-volume-geometry-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/geometries-nurbs-surface-geometry-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/geometries-nurbs-curve-geometry-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/geometries-nurbs-curve-on-surface-geometry-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/geometries-surface-in-nurbs-volume-geometry-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/geometries-brep-surface-h
    resolved_by: tree-sitter
    confidence: exact
  member_of:
  - target: packages/KratosMultiphysics
    resolved_by: tree-sitter
    confidence: exact
---

# Contains

- [SetId1](../../functions/python/add_geometries_to_python/SetId1.md)
- [SetId2](../../functions/python/add_geometries_to_python/SetId2.md)
- [IsIdGeneratedFromString1](../../functions/python/add_geometries_to_python/IsIdGeneratedFromString1.md)
- [IsIdSelfAssigned1](../../functions/python/add_geometries_to_python/IsIdSelfAssigned1.md)
- [Assign](../../functions/python/add_geometries_to_python/Assign.md)
- [Calculate](../../functions/python/add_geometries_to_python/Calculate.md)
- [AddGeometriesToPython](../../functions/python/add_geometries_to_python/AddGeometriesToPython.md)

# Imports

- `"includes/define_python.h"`
- `"geometries/point.h"`
- `"includes/node.h"`
- `"python/containers_interface.h"`
- `"python/add_geometries_to_python.h"`
- `"geometries/geometry.h"`
- `"geometries/line_2d_2.h"`
- `"geometries/line_2d_3.h"`
- `"geometries/triangle_2d_3.h"`
- `"geometries/triangle_2d_6.h"`
- `"geometries/quadrilateral_2d_4.h"`
- `"geometries/quadrilateral_2d_8.h"`
- `"geometries/quadrilateral_2d_9.h"`
- `"geometries/line_3d_2.h"`
- `"geometries/line_3d_3.h"`
- `"geometries/triangle_3d_3.h"`
- `"geometries/triangle_3d_6.h"`
- `"geometries/quadrilateral_3d_4.h"`
- `"geometries/quadrilateral_3d_8.h"`
- `"geometries/quadrilateral_3d_9.h"`
- `"geometries/tetrahedra_3d_4.h"`
- `"geometries/tetrahedra_3d_10.h"`
- `"geometries/pyramid_3d_5.h"`
- `"geometries/pyramid_3d_13.h"`
- `"geometries/prism_3d_6.h"`
- `"geometries/prism_3d_15.h"`
- `"geometries/hexahedra_3d_8.h"`
- `"geometries/hexahedra_3d_20.h"`
- `"geometries/hexahedra_3d_27.h"`
- `"geometries/nurbs_volume_geometry.h"`
- `"geometries/nurbs_surface_geometry.h"`
- `"geometries/nurbs_curve_geometry.h"`
- `"geometries/nurbs_curve_on_surface_geometry.h"`
- `"geometries/surface_in_nurbs_volume_geometry.h"`
- `"geometries/brep_surface.h"`

# Member of

- [KratosMultiphysics](../../packages/KratosMultiphysics.md)