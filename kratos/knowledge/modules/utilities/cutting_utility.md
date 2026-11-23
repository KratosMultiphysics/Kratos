---
type: C++ Module
title: cutting_utility
resource: utilities/cutting_utility.h#L1-L968
generated:
  by: okf-rs/0.7.0
relationships:
  imports:
  - target: external/omp-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/string
    resolved_by: tree-sitter
    confidence: exact
  - target: external/iostream
    resolved_by: tree-sitter
    confidence: exact
  - target: external/cstdlib
    resolved_by: tree-sitter
    confidence: exact
  - target: external/cmath
    resolved_by: tree-sitter
    confidence: exact
  - target: external/algorithm
    resolved_by: tree-sitter
    confidence: exact
  - target: external/boost-numeric-ublas-matrix-hpp
    resolved_by: tree-sitter
    confidence: exact
  - target: external/boost-numeric-ublas-vector-hpp
    resolved_by: tree-sitter
    confidence: exact
  - target: external/boost-numeric-ublas-banded-hpp
    resolved_by: tree-sitter
    confidence: exact
  - target: external/boost-numeric-ublas-matrix-sparse-hpp
    resolved_by: tree-sitter
    confidence: exact
  - target: external/boost-numeric-ublas-triangular-hpp
    resolved_by: tree-sitter
    confidence: exact
  - target: external/boost-numeric-ublas-operation-hpp
    resolved_by: tree-sitter
    confidence: exact
  - target: external/boost-numeric-ublas-lu-hpp
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-define-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-model-part-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-node-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-dof-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-variables-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/containers-array-1d-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/processes-find-nodal-neighbours-process-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/containers-data-value-container-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/includes-mesh-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-math-utils-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-split-tetrahedra-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/utilities-split-triangle-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/geometries-tetrahedra-3d-4-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/geometries-triangle-3d-3-h
    resolved_by: tree-sitter
    confidence: exact
  - target: external/spatial-containers-spatial-containers-h
    resolved_by: tree-sitter
    confidence: exact
  member_of:
  - target: packages/KratosMultiphysics
    resolved_by: tree-sitter
    confidence: exact
---

# Contains

- [CuttingUtility](../../classes/utilities/cutting_utility/CuttingUtility.md)
- [CuttingUtility](../../functions/utilities/cutting_utility/CuttingUtility/CuttingUtility.md)
- [FindSmallestEdge](../../functions/utilities/cutting_utility/CuttingUtility/FindSmallestEdge.md)
- [GenerateCut](../../functions/utilities/cutting_utility/CuttingUtility/GenerateCut.md)
- [AddSkinConditions](../../functions/utilities/cutting_utility/CuttingUtility/AddSkinConditions.md)
- [AddVariablesToCutModelPart](../../functions/utilities/cutting_utility/CuttingUtility/AddVariablesToCutModelPart.md)
- [CSR_Row_Matrix_Mod](../../functions/utilities/cutting_utility/CuttingUtility/CSR_Row_Matrix_Mod.md)
- [FirstLoop](../../functions/utilities/cutting_utility/CuttingUtility/FirstLoop.md)
- [Create_List_Of_New_Nodes_Mod](../../functions/utilities/cutting_utility/CuttingUtility/Create_List_Of_New_Nodes_Mod.md)
- [Calculate_Coordinate_And_Insert_New_Nodes_Mod](../../functions/utilities/cutting_utility/CuttingUtility/Calculate_Coordinate_And_Insert_New_Nodes_Mod.md)
- [GenerateElements](../../functions/utilities/cutting_utility/CuttingUtility/GenerateElements.md)
- [UpdateCutData](../../functions/utilities/cutting_utility/CuttingUtility/UpdateCutData.md)

# Imports

- `<omp.h>`
- `<string>`
- `<iostream>`
- `<cstdlib>`
- `<cmath>`
- `<algorithm>`
- `<boost/numeric/ublas/matrix.hpp>`
- `<boost/numeric/ublas/vector.hpp>`
- `<boost/numeric/ublas/banded.hpp>`
- `<boost/numeric/ublas/matrix_sparse.hpp>`
- `<boost/numeric/ublas/triangular.hpp>`
- `<boost/numeric/ublas/operation.hpp>`
- `<boost/numeric/ublas/lu.hpp>`
- `"includes/define.h"`
- `"includes/model_part.h"`
- `"includes/node.h"`
- `"includes/dof.h"`
- `"includes/variables.h"`
- `"containers/array_1d.h"`
- `"processes/find_nodal_neighbours_process.h"`
- `"containers/data_value_container.h"`
- `"includes/mesh.h"`
- `"utilities/math_utils.h"`
- `"utilities/split_tetrahedra.h"`
- `"utilities/split_triangle.h"`
- `"geometries/tetrahedra_3d_4.h"`
- `"geometries/triangle_3d_3.h"`
- `"spatial_containers/spatial_containers.h"`

# Member of

- [KratosMultiphysics](../../packages/KratosMultiphysics.md)