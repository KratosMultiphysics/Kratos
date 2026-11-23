---
type: C++ Method
title: fill
resource: includes/amatrix_interface.h#L302-L306
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/benchmarks/parallel_utilities_benchmark/RHSElement/CalculateRHS
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/array_1d/array_1d/BOOST_UBLAS_INLINE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/array_1d/array_1d/resize
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/array_1d/array_1d/clear
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/edge_based_data_structure/EdgeBasedDataStructure/CalculateEdgeDataValues
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/bounding_box/BoundingBox/BoundingBox
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/bounding_box/Set
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/line_2d_3/Line2D3/CalculateShapeFunctionsIntegrationPointsLocalGradients
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/line_2d_4/Line2D4/CalculateShapeFunctionsIntegrationPointsLocalGradients
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/line_2d_5/Line2D5/CalculateShapeFunctionsIntegrationPointsLocalGradients
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/line_3d_3/Line3D3/CalculateShapeFunctionsIntegrationPointsLocalGradients
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/nurbs_shape_function_utilities/nurbs_curve_refinement_utilities/NurbsCurveRefinementUtilities/KnotRefinement
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/nurbs_shape_function_utilities/nurbs_curve_refinement_utilities/NurbsCurveRefinementUtilities/DegreeElevation
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/nurbs_shape_function_utilities/nurbs_surface_refinement_utilities/NurbsSurfaceRefinementUtilities/KnotRefinementU
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/nurbs_shape_function_utilities/nurbs_surface_refinement_utilities/NurbsSurfaceRefinementUtilities/KnotRefinementV
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/nurbs_shape_function_utilities/nurbs_surface_refinement_utilities/NurbsSurfaceRefinementUtilities/DegreeElevationU
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/nurbs_shape_function_utilities/nurbs_surface_refinement_utilities/NurbsSurfaceRefinementUtilities/DegreeElevationV
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/sphere_3d_1/Sphere3D1/CalculateShapeFunctionsIntegrationPointsLocalGradients
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/amatrix_interface/Matrix/fill_identity
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/base_64_encoded_output/EncodeTriplet
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/cad_json_output/CadJsonOutput/GetBrepSurfaceParameters
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/cad_json_output/CadJsonOutput/GetBoundaryLoopParameters
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/cartesian_mesh_generator_modeler/CartesianMeshGeneratorModeler/CalculateNormals
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/cartesian_mesh_generator_modeler/CartesianMeshGeneratorModeler/CalculateIsInside
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/tests/test_distributed_sparse_matrices/TestDistributedSparseMatrices/test_matrix_construction
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/reorder_and_optimize_modelpart_process/get
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/assign_average_master_slave_constraints_process/AssignAverageMasterSlaveConstraintsProcess/__AverageMultiPointConstraint
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/p_multigrid/master_slave_constraint_assembler/MasterSlaveConstraintAssembler<TSparse,TDense>/AllocateConstraints
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/p_multigrid/p_multigrid_builder_and_solver/PMultigridBuilderAndSolver<TSparse,TDense>/SetUpDofSet
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/p_multigrid/p_multigrid_builder_and_solver/PMultigridBuilderAndSolver<TSparse,TDense>/ProjectGrid
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/p_multigrid/p_multigrid_utilities/MakePRestrictionOperator
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/spaces/kratos_space/KratosSpace/ResizeData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/spaces/ublas_space/Set
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/spaces/ublas_space/ResizeData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/spaces/ublas_space/SetToZero
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/spatial_containers/bounding_volume_tree/kDOP/Initialize
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/utilities/test_data_type_traits/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/utilities/test_parallel_utilities/RHSElement/CalculateRHS
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_matrix_interface/TestMatrixInterface/test_copy
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_matrix_interface/TestMatrixInterface/test_matrix_vector
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_matrix_interface/TestMatrixInterface/test_matrix_sum
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_matrix_interface/TestMatrixInterface/test_matrix_diff
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_matrix_interface/TestMatrixInterface/test_scalar_prod
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_numpy_export_dense_matrix/TestNumpyExportDenseMatrix/test_numpy_export_dense_matrix_no_copying
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_numpy_export_dense_matrix/TestNumpyExportDenseMatrix/test_numpy_export_complex_dense_matrix_no_copying
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_numpy_export_dense_matrix/TestNumpyExportDenseMatrix/test_numpy_export_dense_matrix_copying
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_vector_interface/TestVectorInterface/test_sum
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/entities_utilities/EntitityIdentifier<TEntity>/EntitityIdentifier
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/search_utilities/SynchronizePointsWithRanks
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/sparse_matrix_multiplication_utility/TLS/TLS
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/sparse_matrix_multiplication_utility/MatrixAdd
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void fill(data_type const& value)`

# Called by

- [CalculateRHS](../../../../functions/benchmarks/parallel_utilities_benchmark/RHSElement/CalculateRHS.md)
- [BOOST_UBLAS_INLINE](../../../../functions/containers/array_1d/array_1d/BOOST_UBLAS_INLINE.md)
- [resize](../../../../functions/containers/array_1d/array_1d/resize.md)
- [clear](../../../../functions/containers/array_1d/array_1d/clear.md)
- [CalculateEdgeDataValues](../../../../functions/containers/edge_based_data_structure/EdgeBasedDataStructure/CalculateEdgeDataValues.md)
- [BoundingBox](../../../../functions/geometries/bounding_box/BoundingBox/BoundingBox.md)
- [Set](../../../../functions/geometries/bounding_box/Set.md)
- [CalculateShapeFunctionsIntegrationPointsLocalGradients](../../../../functions/geometries/line_2d_3/Line2D3/CalculateShapeFunctionsIntegrationPointsLocalGradients.md)
- [CalculateShapeFunctionsIntegrationPointsLocalGradients](../../../../functions/geometries/line_2d_4/Line2D4/CalculateShapeFunctionsIntegrationPointsLocalGradients.md)
- [CalculateShapeFunctionsIntegrationPointsLocalGradients](../../../../functions/geometries/line_2d_5/Line2D5/CalculateShapeFunctionsIntegrationPointsLocalGradients.md)
- [CalculateShapeFunctionsIntegrationPointsLocalGradients](../../../../functions/geometries/line_3d_3/Line3D3/CalculateShapeFunctionsIntegrationPointsLocalGradients.md)
- [KnotRefinement](../../../../functions/geometries/nurbs_shape_function_utilities/nurbs_curve_refinement_utilities/NurbsCurveRefinementUtilities/KnotRefinement.md)
- [DegreeElevation](../../../../functions/geometries/nurbs_shape_function_utilities/nurbs_curve_refinement_utilities/NurbsCurveRefinementUtilities/DegreeElevation.md)
- [KnotRefinementU](../../../../functions/geometries/nurbs_shape_function_utilities/nurbs_surface_refinement_utilities/NurbsSurfaceRefinementUtilities/KnotRefinementU.md)
- [KnotRefinementV](../../../../functions/geometries/nurbs_shape_function_utilities/nurbs_surface_refinement_utilities/NurbsSurfaceRefinementUtilities/KnotRefinementV.md)
- [DegreeElevationU](../../../../functions/geometries/nurbs_shape_function_utilities/nurbs_surface_refinement_utilities/NurbsSurfaceRefinementUtilities/DegreeElevationU.md)
- [DegreeElevationV](../../../../functions/geometries/nurbs_shape_function_utilities/nurbs_surface_refinement_utilities/NurbsSurfaceRefinementUtilities/DegreeElevationV.md)
- [CalculateShapeFunctionsIntegrationPointsLocalGradients](../../../../functions/geometries/sphere_3d_1/Sphere3D1/CalculateShapeFunctionsIntegrationPointsLocalGradients.md)
- [fill_identity](../../../../functions/includes/amatrix_interface/Matrix/fill_identity.md)
- [EncodeTriplet](../../../../functions/input_output/base_64_encoded_output/EncodeTriplet.md)
- [GetBrepSurfaceParameters](../../../../functions/input_output/cad_json_output/CadJsonOutput/GetBrepSurfaceParameters.md)
- [GetBoundaryLoopParameters](../../../../functions/input_output/cad_json_output/CadJsonOutput/GetBoundaryLoopParameters.md)
- [CalculateNormals](../../../../functions/modeler/cartesian_mesh_generator_modeler/CartesianMeshGeneratorModeler/CalculateNormals.md)
- [CalculateIsInside](../../../../functions/modeler/cartesian_mesh_generator_modeler/CartesianMeshGeneratorModeler/CalculateIsInside.md)
- [test_matrix_construction](../../../../functions/mpi/tests/test_distributed_sparse_matrices/TestDistributedSparseMatrices/test_matrix_construction.md)
- [get](../../../../functions/processes/reorder_and_optimize_modelpart_process/get.md)
- [__AverageMultiPointConstraint](../../../../functions/python_scripts/assign_average_master_slave_constraints_process/AssignAverageMasterSlaveConstraintsProcess/__AverageMultiPointConstraint.md)
- [AllocateConstraints](../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/master_slave_constraint_assembler/MasterSlaveConstraintAssembler<TSparse,TDense>/AllocateConstraints.md)
- [SetUpDofSet](../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/p_multigrid_builder_and_solver/PMultigridBuilderAndSolver<TSparse,TDense>/SetUpDofSet.md)
- [ProjectGrid](../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/p_multigrid_builder_and_solver/PMultigridBuilderAndSolver<TSparse,TDense>/ProjectGrid.md)
- [MakePRestrictionOperator](../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/p_multigrid_utilities/MakePRestrictionOperator.md)
- [ResizeData](../../../../functions/spaces/kratos_space/KratosSpace/ResizeData.md)
- [Set](../../../../functions/spaces/ublas_space/Set.md)
- [ResizeData](../../../../functions/spaces/ublas_space/ResizeData.md)
- [SetToZero](../../../../functions/spaces/ublas_space/SetToZero.md)
- [Initialize](../../../../functions/spatial_containers/bounding_volume_tree/kDOP/Initialize.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/utilities/test_data_type_traits/KRATOS_TEST_CASE_IN_SUITE.md)
- [CalculateRHS](../../../../functions/tests/cpp_tests/utilities/test_parallel_utilities/RHSElement/CalculateRHS.md)
- [test_copy](../../../../functions/tests/test_matrix_interface/TestMatrixInterface/test_copy.md)
- [test_matrix_vector](../../../../functions/tests/test_matrix_interface/TestMatrixInterface/test_matrix_vector.md)
- [test_matrix_sum](../../../../functions/tests/test_matrix_interface/TestMatrixInterface/test_matrix_sum.md)
- [test_matrix_diff](../../../../functions/tests/test_matrix_interface/TestMatrixInterface/test_matrix_diff.md)
- [test_scalar_prod](../../../../functions/tests/test_matrix_interface/TestMatrixInterface/test_scalar_prod.md)
- [test_numpy_export_dense_matrix_no_copying](../../../../functions/tests/test_numpy_export_dense_matrix/TestNumpyExportDenseMatrix/test_numpy_export_dense_matrix_no_copying.md)
- [test_numpy_export_complex_dense_matrix_no_copying](../../../../functions/tests/test_numpy_export_dense_matrix/TestNumpyExportDenseMatrix/test_numpy_export_complex_dense_matrix_no_copying.md)
- [test_numpy_export_dense_matrix_copying](../../../../functions/tests/test_numpy_export_dense_matrix/TestNumpyExportDenseMatrix/test_numpy_export_dense_matrix_copying.md)
- [test_sum](../../../../functions/tests/test_vector_interface/TestVectorInterface/test_sum.md)
- [EntitityIdentifier](../../../../functions/utilities/entities_utilities/EntitityIdentifier<TEntity>/EntitityIdentifier.md)
- [SynchronizePointsWithRanks](../../../../functions/utilities/search_utilities/SynchronizePointsWithRanks.md)
- [TLS](../../../../functions/utilities/sparse_matrix_multiplication_utility/TLS/TLS.md)
- [MatrixAdd](../../../../functions/utilities/sparse_matrix_multiplication_utility/MatrixAdd.md)