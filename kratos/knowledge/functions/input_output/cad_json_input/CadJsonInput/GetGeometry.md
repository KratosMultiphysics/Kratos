---
type: C++ Method
title: GetGeometry
resource: input_output/cad_json_input.h#L931-L941
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/GetInt
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetString
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/conditions/mesh_condition/MeshCondition/Create
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/conditions/mesh_condition/MeshCondition/Clone
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/conditions/mesh_condition/MeshCondition/CalculateOnIntegrationPoints
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/edge_based_data_structure/EdgeBasedDataStructure/CalculateEdgeDataValues
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/elements/distance_calculation_element_simplex/DistanceCalculationElementSimplex/Create
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/elements/distance_calculation_element_simplex/DistanceCalculationElementSimplex/CalculateLocalSystem
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/elements/distance_calculation_element_simplex/DistanceCalculationElementSimplex/EquationIdVector
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/elements/distance_calculation_element_simplex/DistanceCalculationElementSimplex/GetDofList
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/elements/distance_calculation_element_simplex/DistanceCalculationElementSimplex/Check
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/elements/edge_based_gradient_recovery_element/EdgeBasedGradientRecoveryElement<TDim>/CalculateLocalSystem
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/elements/edge_based_gradient_recovery_element/EdgeBasedGradientRecoveryElement<TDim>/Check
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/elements/edge_based_gradient_recovery_element/EdgeBasedGradientRecoveryElement<TDim>/CalculateLeftHandSide
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/elements/edge_based_gradient_recovery_element/EdgeBasedGradientRecoveryElement<TDim>/CalculateRightHandSide
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/elements/edge_based_gradient_recovery_element/EdgeBasedGradientRecoveryElement<2>/EquationIdVector
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/elements/edge_based_gradient_recovery_element/EdgeBasedGradientRecoveryElement<3>/EquationIdVector
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/elements/edge_based_gradient_recovery_element/EdgeBasedGradientRecoveryElement<2>/GetDofList
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/elements/edge_based_gradient_recovery_element/EdgeBasedGradientRecoveryElement<3>/GetDofList
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/elements/edge_based_gradient_recovery_element/Create
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/elements/embedded_nodal_variable_calculation_element_simplex/EmbeddedNodalVariableCalculationElementSimplex<double>/Check
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/elements/embedded_nodal_variable_calculation_element_simplex/EmbeddedNodalVariableCalculationElementSimplex<array_1d<double,3>>/Check
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/elements/embedded_nodal_variable_calculation_element_simplex/EmbeddedNodalVariableCalculationElementSimplex<double>/CalculateRightHandSide
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/elements/embedded_nodal_variable_calculation_element_simplex/EmbeddedNodalVariableCalculationElementSimplex<array_1d<double, 3>>/CalculateRightHandSide
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/elements/embedded_nodal_variable_calculation_element_simplex/EmbeddedNodalVariableCalculationElementSimplex<double>/EquationIdVector
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/elements/embedded_nodal_variable_calculation_element_simplex/EmbeddedNodalVariableCalculationElementSimplex<array_1d<double, 3>>/EquationIdVector
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/elements/embedded_nodal_variable_calculation_element_simplex/EmbeddedNodalVariableCalculationElementSimplex<double>/GetDofList
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/elements/embedded_nodal_variable_calculation_element_simplex/EmbeddedNodalVariableCalculationElementSimplex<array_1d<double, 3>>/GetDofList
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/elements/embedded_nodal_variable_calculation_element_simplex/Create
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/elements/levelset_convection_element_simplex/LevelSetConvectionElementSimplex/Create
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/elements/levelset_convection_element_simplex/LevelSetConvectionElementSimplex/Initialize
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/elements/levelset_convection_element_simplex/LevelSetConvectionElementSimplex/CalculateLocalSystem
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/elements/levelset_convection_element_simplex/LevelSetConvectionElementSimplex/EquationIdVector
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/elements/levelset_convection_element_simplex/LevelSetConvectionElementSimplex/GetDofList
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/elements/levelset_convection_element_simplex_algebraic_stabilization/LevelSetConvectionElementSimplexAlgebraicStabilization<TDim, TNumNodes>/Create
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/elements/levelset_convection_element_simplex_algebraic_stabilization/LevelSetConvectionElementSimplexAlgebraicStabilization<TDim, TNumNodes>/CalculateLocalSystem
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/elements/mesh_element/MeshElement/Create
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/elements/mesh_element/MeshElement/Clone
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/elements/mesh_element/MeshElement/CalculateOnIntegrationPoints
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/future/tests/cpp_tests/test_utilities/solving_strategies_test_utilities/TestLaplacianElement1D/CalculateLeftHandSide
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/future/tests/cpp_tests/test_utilities/solving_strategies_test_utilities/TestLaplacianElement1D/CalculateRightHandSide
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/future/tests/cpp_tests/test_utilities/solving_strategies_test_utilities/TestLaplacianElement1D/EquationIdVector
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/future/tests/cpp_tests/test_utilities/solving_strategies_test_utilities/TestLaplacianElement1D/GetDofList
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/future/tests/cpp_tests/test_utilities/solving_strategies_test_utilities/TestElasticityElement2D/CalculateLeftHandSide
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/future/tests/cpp_tests/test_utilities/solving_strategies_test_utilities/TestElasticityElement2D/CalculateRightHandSide
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/future/tests/cpp_tests/test_utilities/solving_strategies_test_utilities/TestElasticityElement2D/EquationIdVector
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/future/tests/cpp_tests/test_utilities/solving_strategies_test_utilities/TestElasticityElement2D/GetDofList
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/condition/Condition/Create
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/condition/Condition/Clone
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/condition/Condition/Check
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/element/Element/Create
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/element/Element/Clone
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/element/Element/Check
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/geometrical_object/SetData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/geometrical_object/HasSameGeometryType
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/gid_gauss_point_container/GidGaussPointsContainer/AddElement
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/gid_gauss_point_container/GidGaussPointsContainer/AddCondition
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/gid_io/WriteSphereMesh
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/gid_io/WriteClusterMesh
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/gid_mesh_container/GidMeshContainer/AddElement
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/gid_mesh_container/GidMeshContainer/AddCondition
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/gid_mesh_container/GidMeshContainer/WriteMesh
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/mesh/Mesh/WorkingSpaceDimension
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/cad_json_input/CadJsonInput/ReadBrepEdgeBrepCurveOnSurface
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/cad_json_input/CadJsonInput/ReadCouplingGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/cad_json_input/CadJsonInput/ReadPointsOnGeometries
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/ensight_output/EnSightOutput/WriteGeometricalGaussVariableToFile
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/ensight_output/EnSightOutput/CollectPartData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/ensight_output/if
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/ensight_output/EnSightOutput/GetGeometryConnectivity
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/ensight_output/EnSightOutput/GetAverageIntegrationValue
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/obj_io/ObjIO/WriteModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/stl_io/StlIO/WriteEntityBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/stl_io/StlIO/WriteEntityBlockMPI
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/unv_output/UnvOutput/WriteElements
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/vtk_output/VtkOutput/DetermineVtkCellListSize
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/vtk_output/VtkOutput/WriteCellType
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/vtk_output/VtkOutput/WriteIntegrationScalarContainerVariable
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/vtk_output/VtkOutput/WriteIntegrationVectorContainerVariable
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/vtk_output/VtkOutput/WriteModelPartWithoutNodesToFile
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/vtu_output/GetGeometryTypes
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/vtu_output/GetOffsets
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/vtu_output/GetConnectivities
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/vtu_output/GetNodesContainer
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/vtu_output/VtuOutput/WriteIntegrationPointData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/cartesian_mesh_generator_modeler/CartesianMeshGeneratorModeler/GenerateMesh
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/cartesian_mesh_generator_modeler/CartesianMeshGeneratorModeler/CalculateNormals
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/cartesian_mesh_generator_modeler/CartesianMeshGeneratorModeler/CalculateBoundaryIntersections
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/cartesian_mesh_generator_modeler/CartesianMeshGeneratorModeler/CalculateBoundingBox
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/clean_up_problematic_triangles_modeler/CleanUpProblematicTrianglesModeler/CleanUpProblematicGeometries
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/clean_up_problematic_triangles_modeler/CleanUpProblematicTrianglesModeler/ComputeNullAreaTriangles
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/coloring/color_and_intersect_with_cells_in_touch/ColorAndIntersectWithCellsInTouch/Apply
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/coloring/color_cell_faces/ColorCellFaces/Apply
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/coloring/color_cells_in_touch/ColorCellsInTouch/Apply
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/coloring/color_cells_with_inside_center/ColorCellsWithInsideCenter/Apply
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/coloring/color_connected_cells_in_touch/ColorConnectedCellsInTouch/Apply
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/duplicate_mesh_modeler/DuplicateMeshModeler/GenerateMesh
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/edge_swapping_2d_modeler/EdgeSwapping2DModeler/WriteMesh
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/edge_swapping_2d_modeler/EdgeSwapping2DModeler/CollapseNodes
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/edge_swapping_2d_modeler/EdgeSwapping2DModeler/SetNodalCollapsingData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/edge_swapping_2d_modeler/EdgeSwapping2DModeler/Swap
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/edge_swapping_2d_modeler/EdgeSwapping2DModeler/ElementQuality
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/edge_swapping_2d_modeler/EdgeSwapping2DModeler/FindNodalNeighbours
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/edge_swapping_2d_modeler/EdgeSwapping2DModeler/SetSwappingData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/edge_swapping_2d_modeler/EdgeSwapping2DModeler/FindNeighbourElement
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/edge_swapping_2d_modeler/EdgeSwapping2DModeler/IsInElementSphere
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/entity_generation/generate_stl_intersection_with_cells/GenerateStlIntersectionWithCells/GenerateIntersection
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/mesh_suite_modeler/MeshSuiteModeler/SetMesh
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/operation/compute_surrogate_boundary_data/ComputeSurrogateBoundaryData/Execute
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/operation/find_contacts_in_skin_model_part/FindContactsInSkinModelPart/FindConditionContact
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/structured_mesh_refinement_modeler/StructuredMeshRefinementModeler/GenerateMesh
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/structured_mesh_refinement_modeler/GenerateNodes
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/structured_mesh_refinement_modeler/GenerateNode
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/structured_mesh_refinement_modeler/StructuredMeshRefinementModeler/GenerateElements
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/structured_mesh_refinement_modeler/StructuredMeshRefinementModeler/GenerateConditions
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/future/tests/test_mpi_vtu_output/TestMPIVtuOutput/__InitializeModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/tests/test_mpi_serializer/executeComputeArea_Task
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/utilities/gather_modelpart_utility/GatherModelPartUtility/GatherModelPartUtility
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/utilities/mpi_normal_calculation_utilities/MPINormalCalculationUtils/OrientFaces
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/utilities/mpi_normal_calculation_utilities/MPINormalCalculationUtils/IdentifyFaces
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/utilities/mpi_normal_calculation_utilities/MPINormalCalculationUtils/InitializeNormalData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/utilities/mpi_normal_calculation_utilities/MPINormalCalculationUtils/UpdateNodeNormals
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/apply_periodic_boundary_condition_process/ApplyPeriodicConditionProcess/ApplyConstraintsForPeriodicConditions
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/apply_ray_casting_process/ApplyRayCastingProcess<TDim>/GetCellIntersections
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/assign_scalar_field_to_entities_process/AssignScalarFieldToEntitiesProcess<Condition>/CallFunction
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/assign_scalar_field_to_entities_process/AssignScalarFieldToEntitiesProcess<Element>/CallFunction
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/assign_scalar_field_to_entities_process/AssignScalarFieldToEntitiesProcess<Condition>/CallFunctionComponents
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/assign_scalar_field_to_entities_process/AssignScalarFieldToEntitiesProcess<Element>/CallFunctionComponents
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/assign_scalar_field_to_entities_process/AssignScalarFieldToEntitiesProcess<Condition>/CallFunctionLocalSystem
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/assign_scalar_field_to_entities_process/AssignScalarFieldToEntitiesProcess<Element>/CallFunctionLocalSystem
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/assign_scalar_field_to_entities_process/AssignScalarFieldToEntitiesProcess<Condition>/CallFunctionLocalSystemComponents
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/assign_scalar_field_to_entities_process/AssignScalarFieldToEntitiesProcess<Element>/CallFunctionLocalSystemComponents
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/assign_scalar_field_to_entities_process/AssignScalarFieldToEntitiesProcess<Condition>/AssignTimeDependentValue
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/assign_scalar_field_to_entities_process/AssignScalarFieldToEntitiesProcess<Element>/AssignTimeDependentValue
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/assign_scalar_input_to_entities_process/AssignScalarInputToEntitiesProcess<Condition, AssignScalarInputToEntitiesProcessSettings/SaveAsNonHistoricalVariable>/GetCoordinatesEntity
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/assign_scalar_input_to_entities_process/AssignScalarInputToEntitiesProcess<Element, AssignScalarInputToEntitiesProcessSettings/SaveAsNonHistoricalVariable>/GetCoordinatesEntity
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_discontinuous_distance_to_skin_process/CalculateDiscontinuousDistanceToSkinProcess<TDim>/CalculateElementalDistances
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_discontinuous_distance_to_skin_process/CalculateDiscontinuousDistanceToSkinProcess<TDim>/CalculateElementalAndEdgeDistances
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_discontinuous_distance_to_skin_process/CalculateDiscontinuousDistanceToSkinProcess<TDim>/ComputeEdgesIntersections
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_discontinuous_distance_to_skin_process/CalculateDiscontinuousDistanceToSkinProcess<TDim>/ComputeElementalDistancesFromPlaneApproximation
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_discontinuous_distance_to_skin_process/CalculateDiscontinuousDistanceToSkinProcess<TDim>/ComputeIntersectionPlaneElementalDistances
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_discontinuous_distance_to_skin_process/CalculateDiscontinuousDistanceToSkinProcess<TDim>/ComputePlaneApproximation
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_discontinuous_distance_to_skin_process/CalculateDiscontinuousDistanceToSkinProcess<TDim>/CorrectDistanceOrientation
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_discontinuous_distance_to_skin_process/CalculateDiscontinuousDistanceToSkinProcess<TDim>/ComputeElementalDistancesFromEdgeRatios
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_discontinuous_distance_to_skin_process/CalculateDiscontinuousDistanceToSkinProcess<TDim>/CheckAndCorrectEdgeDistances
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_discontinuous_distance_to_skin_process/CalculateEmbeddedVariableFromSkinSpecialization
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_distance_to_path_process/CalculateDistanceToPathProcess<THistorical>/Execute
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_distance_to_skin_process/CalculateDistanceToSkinProcess<TDim>/CalculateElementalDistances
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_distance_to_skin_process/CalculateDistanceToSkinProcess<TDim>/CalculateDistanceToNode
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_distance_to_skin_process/CalculateDistanceToSkinProcess<TDim>/CalculateNodalDistances
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_embedded_nodal_variable_from_skin_process/CalculateEmbeddedNodalVariableFromSkinProcess/CalculateEmbeddedNodalVariableFromSkinProcess
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_embedded_nodal_variable_from_skin_process/CalculateEmbeddedNodalVariableFromSkinProcess/AddIntersectedElementsModelPartElements
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_embedded_signed_distance_to_3d_skin_process/CalculateEmbeddedSignedDistanceTo3DSkinProcess/CalculateEmbeddedVelocity
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_embedded_signed_distance_to_3d_skin_process/CalculateEmbeddedSignedDistanceTo3DSkinProcess/SetMaximumAndMinimumDistanceValues
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_nodal_area_process/CalculateNodalAreaProcess<THistorical>/Execute
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_signed_distance_to_3d_condition_skin_process/DistanceSpatialContainersConditionConfigure/CalculateBoundingBox
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_signed_distance_to_3d_condition_skin_process/DistanceSpatialContainersConditionConfigure/GetBoundingBox
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_signed_distance_to_3d_condition_skin_process/DistanceSpatialContainersConditionConfigure/Intersection
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_signed_distance_to_3d_condition_skin_process/DistanceSpatialContainersConditionConfigure/IntersectionBox
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/IdentifyIntersectionNodes
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/CalcSignedDistancesToOneIntNode
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/CalcSignedDistancesToTwoIntNodes
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/CalcSignedDistancesToApproxTriangle
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/AssignDistancesToElements
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/CalculateNotEmptyLeavesDistance
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/GetCellIntersections
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_signed_distance_to_3d_skin_process/DistanceSpatialContainersConfigure/CalculateBoundingBox
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_signed_distance_to_3d_skin_process/DistanceSpatialContainersConfigure/GetBoundingBox
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_signed_distance_to_3d_skin_process/DistanceSpatialContainersConfigure/Intersection
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_signed_distance_to_3d_skin_process/DistanceSpatialContainersConfigure/IntersectionBox
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/MappingPressureToStructure
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/AveragePressureToNode
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/IdentifyIntersectionNodes
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/CalcSignedDistancesToOneIntNode
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/CalcSignedDistancesToTwoIntNodes
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/CalcSignedDistancesToThreeIntNodes
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/CalcSignedDistancesToMoreThanThreeIntNodes
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/AssignMinimalNodalDistance
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/GenerateSkinModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/CalculateNotEmptyLeavesDistance
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/GetCellIntersections
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/check_same_modelpart_using_skin_distance_process/CheckSameModelPartUsingSkinDistanceProcess<TDim>/Execute
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/check_skin_process/CheckSkinProcess/Execute
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/compute_nodal_gradient_process/ComputeNodalGradientProcess<TOutputHistorical>/ComputeElementalContributionsAndVolume
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/compute_nodal_normal_divergence_process/ComputeNodalNormalDivergenceProcess<THistorical>/Execute
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/edge_based_gradient_recovery_process/EdgeBasedGradientRecoveryProcess/CheckAndAssignSettings
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/fast_transfer_between_model_parts_process/FastTransferBetweenModelPartsProcess/ReplicateWithoutFlags
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/fast_transfer_between_model_parts_process/FastTransferBetweenModelPartsProcess/ReplicateWithFlags
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/find_conditions_neighbours_process/FindConditionsNeighboursProcess/Execute
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/find_conditions_neighbours_process/FindConditionsNeighboursProcess/CheckForNeighbourFaces
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/find_conditions_neighbours_process/FindConditionsNeighboursProcess/ComputeDimension
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/find_global_nodal_entity_neighbours_process/FindGlobalNodalEntityNeighboursProcess<TContainerType>/Execute
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/find_global_nodal_neighbours_for_entities_process/FindNodalNeighboursForEntitiesProcess<TContainerType>/Execute
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/find_intersected_geometrical_objects_process/FindIntersectedGeometricalObjectsProcess/WorkingSpaceDimension
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/find_intersected_geometrical_objects_process/FindIntersectedGeometricalObjectsProcess/MarkIfIntersected
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/find_intersected_geometrical_objects_process/FindIntersectedGeometricalObjectsProcess/FindIntersectedSkinObjects
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/find_intersected_geometrical_objects_process/DistanceSpatialContainersConfigure/CalculateBoundingBox
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/find_intersected_geometrical_objects_process/DistanceSpatialContainersConfigure/GetBoundingBox
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/find_intersected_geometrical_objects_process/DistanceSpatialContainersConfigure/Intersection
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/find_intersected_geometrical_objects_process/DistanceSpatialContainersConfigure/IntersectionBox
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/find_nodal_h_process/FindNodalHProcess<THistorical>/Execute
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/generic_find_elements_neighbours_process/GenericFindElementalNeighboursProcess/Execute
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/integration_values_extrapolation_to_nodes_process/IntegrationValuesExtrapolationToNodesProcess/ExecuteFinalizeSolutionStep
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/integration_values_extrapolation_to_nodes_process/IntegrationValuesExtrapolationToNodesProcess/InitializeMaps
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/levelset_convection_process/LevelSetConvectionProcess/EvaluateNumberOfSubsteps
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/levelset_convection_process/LevelSetConvectionProcess/EvaluateLimiter
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/levelset_convection_process/LevelSetConvectionProcess/InitializeConvectionStrategy
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/measure_mesh_quality_process/MeasureMeshQualityProcess/Execute
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/mesh_local_smoothing_process/MeshLocalSmoothingProcess/SelectLowQualityElementNodes
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/mesh_node_collapsing_process/MeshNodeCollapsingProcess/CollapseNode
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/mesh_node_collapsing_process/MeshNodeCollapsingProcess/CalculateMinQualityOfNeighbourElements
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/mesh_node_collapsing_process/MeshNodeCollapsingProcess/ElementHas
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/mesh_node_collapsing_process/MeshNodeCollapsingProcess/SwapElementNode
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/parallel_distance_calculation_process/ParallelDistanceCalculationProcess<TDim>/FindMaximumEdgeSize
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/parallel_distance_calculation_process/ParallelDistanceCalculationProcess<TDim>/SetDistancesOnDividedElements
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/parallel_distance_calculation_process/ParallelDistanceCalculationProcess<TDim>/ExtendDistancesByLayer
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/reorder_and_optimize_modelpart_process/ReorderAndOptimizeModelPartProcess/Execute
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/reorder_and_optimize_modelpart_process/ReorderAndOptimizeModelPartProcess/OptimizeOrdering
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/reorder_and_optimize_modelpart_process/ReorderAndOptimizeModelPartProcess/ReorderElements
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/set_initial_state_process/SetInitialStateProcess<TDim>/ExecuteInitializeSolutionStep
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/simple_mortar_mapper_process/SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>/GetReferenceArea
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/simple_mortar_mapper_process/PerformMortarOperations
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/simple_mortar_mapper_process/ClearIndexes
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/simple_mortar_mapper_process/FillDatabase
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/simple_mortar_mapper_wrapper_process/SimpleMortarMapperProcessWrapper/SimpleMortarMapperProcessWrapper
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/skin_detection_process/SkinDetectionProcess<TDim>/GenerateFaceMaps
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/skin_detection_process/SkinDetectionProcess<TDim>/SetUpAdditionalSubModelParts
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/split_internal_interfaces_process/SplitInternalInterfacesProcess/ExecuteInitialize
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/split_internal_interfaces_process/SplitInternalInterfacesProcess/SplitBoundary
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/tetrahedra_mesh_edge_swapping_process/TetrahedraMeshEdgeSwappingProcess/Execute
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/tetrahedra_mesh_quality_weighted_smoothing_process/TetrahedraMeshQualityWeightedSmoothingProcess/FindOptimumPositionsAndWeights
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/tetrahedra_mesh_worst_element_smoothing_process/TetrahedraMeshWorstElementSmoothingProcess/FindOptimumPositionsAndWeights
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/tetrahedral_mesh_orientation_check/TetrahedralMeshOrientationCheck/Execute
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/tetrahedral_mesh_orientation_check/TetrahedralMeshOrientationCheck/SwapAll
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/tetrahedral_mesh_orientation_check/TetrahedralMeshOrientationCheck/SwapNegativeElements
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/variational_distance_calculation_process/VariationalDistanceCalculationProcess/Execute
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/variational_distance_calculation_process/VariationalDistanceCalculationProcess/ValidateInput
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/variational_distance_calculation_process/VariationalDistanceCalculationProcess/ReGenerateDistanceModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python/add_geometrical_utilities_to_python/VectorizedFindHelper
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python/add_mesh_to_python/GetNodeFromElement
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python/add_mesh_to_python/GetNodesFromElement
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python/add_mesh_to_python/GetNodeFromCondition
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python/add_mesh_to_python/GetNodesFromCondition
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python/add_mesh_to_python/GetIntegrationPointsFromElement
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python/add_mesh_to_python/SetValuesOnIntegrationPoints
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python/add_mesh_to_python/SetValuesOnIntegrationPointsVector
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python/add_mesh_to_python/SetValuesOnIntegrationPointsConstitutiveLaw
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/p_multigrid/p_multigrid_utilities/MakePRestrictionOperator
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_componentwise/ResidualBasedEliminationBuilderAndSolverComponentwise/Build
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/ResidualBasedAdjointBossakScheme/CalculatePreviousTimeStepContributions
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/ResidualBasedAdjointBossakScheme/CalculateNodeNeighbourCount
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/UpdateEntityTimeSchemeContributions
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/UpdateEntityAuxiliaryVariableContributions
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/schemes/residualbased_incrementalupdate_static_scheme_slip/ResidualBasedIncrementalUpdateStaticSchemeSlip/CalculateSystemContributions
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/schemes/residualbased_incrementalupdate_static_scheme_slip/ResidualBasedIncrementalUpdateStaticSchemeSlip/CalculateRHSContribution
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/schemes/residualbased_incrementalupdate_static_scheme_slip/ResidualBasedIncrementalUpdateStaticSchemeSlip/CalculateLHSContribution
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/schemes/sensitivity_builder_scheme/CalculateLocalSensitivityAndGlobalPointersVector
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/WriteElements
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/WriteConditions
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/FillNodalConnectivitiesFromElementBlockInList
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/FillNodalConnectivitiesFromConditionBlockInList
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadElementsBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadConditionsBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadElementsConnectivitiesBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadConditionsConnectivitiesBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/FillNodalConnectivitiesFromElementBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/FillNodalConnectivitiesFromConditionBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/DivideElementsBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/DivideConditionsBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/periodic_condition/PeriodicCondition/Create
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/periodic_condition/PeriodicCondition/GetDofList
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/spatial_containers/bounding_volume_tree/LineRegressionVolumePartitioner/Partition
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/spatial_containers/bounding_volume_tree/SimpleBoundingVolumePartitioner/Partition
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/spatial_containers/bounding_volume_tree/BoundingVolumeTree/UpdateTree
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/spatial_containers/bounding_volume_tree/BoundingVolumeTree/AddGeometryFromConditions
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/spatial_containers/dem_search/RadiusPoint/Initialize
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/spatial_containers/geometrical_objects_bins/GeometricalObjectsBins/SearchInRadius
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/spatial_containers/geometrical_objects_bins/GeometricalObjectsBins/SearchNearestInCell
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/spatial_containers/geometrical_objects_bins/GeometricalObjectsBins/SearchIsInsideInCell
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/spatial_containers/geometrical_objects_bins/GeometricalObjectsBins
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/spatial_containers/geometrical_objects_bins/AddObjectsToCells
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/spatial_containers/octree_binary/OctreeBinary/IsIntersected
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/spatial_containers/point_object/PointObject<TObject>/UpdatePoint
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/spatial_containers/stl_bvh_tree/StlBvhTree
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/templates/classes/condition_template/EquationIdVector
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/templates/classes/condition_template/GetDofList
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/templates/classes/condition_template/Check
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/templates/classes/element_template/EquationIdVector
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/templates/classes/element_template/GetDofList
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/templates/classes/element_template/Check
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tensor_adaptors/connectivity_ids_tensor_adaptor/ConnectivityIdsTensorAdaptor/GetShape
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tensor_adaptors/connectivity_ids_tensor_adaptor/ConnectivityIdsTensorAdaptor/CheckContainer
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tensor_adaptors/connectivity_ids_tensor_adaptor/ConnectivityIdsTensorAdaptor/CollectIds
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tensor_adaptors/geometries_tensor_adaptor/GeometriesTensorAdaptor/GetIntegrationMethod
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tensor_adaptors/geometries_tensor_adaptor/GeometriesTensorAdaptor/GetShape
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tensor_adaptors/geometries_tensor_adaptor/GeometriesTensorAdaptor/CollectIntegrationWeights
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tensor_adaptors/geometries_tensor_adaptor/GeometriesTensorAdaptor/CollectShapeFunctions
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tensor_adaptors/geometries_tensor_adaptor/GeometriesTensorAdaptor/CollectShapeFunctionsDerivatives
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tensor_adaptors/geometries_tensor_adaptor/GeometriesTensorAdaptor/CollectNormals
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tensor_adaptors/geometries_tensor_adaptor/GeometriesTensorAdaptor/CollectJacobians
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tensor_adaptors/geometry_metrics_tensor_adaptor/FillDomainSize
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/containers/test_model_part_geometry_container/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/elements/test_edge_based_gradient_recovery_element/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/elements/test_levelset_convection_element_simplex/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/external_libraries/test_delaunator/validatewithtolerance
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/geometries/test_coupling_geometry/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/geometries/test_nurbs_curve/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/geometries/test_nurbs_surface/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/geometries/test_nurbs_volume/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/modeler/test_voxel_mesh_modeler/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/processes/test_calculate_discontinuous_distance_to_skin_process/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/processes/test_coarsening_process/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/processes/test_find_nodal_neighbours_process/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/processes/test_mortar_mapper_process/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/processes/test_parallel_distance_calculation_process/SetUpDistanceField
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/processes/test_replace_elements_and_condition_process/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/processes/test_set_initial_state_process/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/processes/test_structured_mesh_generator_process/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/sources/test_condition/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/PrimalElement/EquationIdVector
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/PrimalElement/GetDofList
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/PrimalElement/GetValuesVector
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/PrimalElement/GetFirstDerivativesVector
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/PrimalElement/GetSecondDerivativesVector
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/PrimalElement/CalculateLeftHandSide
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/PrimalElement/CalculateRightHandSide
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/PrimalCondition/EquationIdVector
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/PrimalCondition/GetDofList
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/PrimalCondition/GetValuesVector
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/PrimalCondition/GetFirstDerivativesVector
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/PrimalCondition/GetSecondDerivativesVector
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/PrimalCondition/CalculateLeftHandSide
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/PrimalCondition/CalculateRightHandSide
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/ThisExtensions/GetFirstDerivativesVector
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/ThisExtensions/GetSecondDerivativesVector
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/ThisExtensions/GetAuxiliaryVector
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/AdjointElement/EquationIdVector
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/AdjointElement/GetDofList
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/AdjointElement/GetValuesVector
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/AdjointElement/GetFirstDerivativesVector
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/AdjointElement/GetSecondDerivativesVector
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/AdjointElement/CalculateSensitivityMatrix
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/AdjointCondition/EquationIdVector
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/AdjointCondition/GetDofList
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/AdjointCondition/GetValuesVector
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/AdjointCondition/GetFirstDerivativesVector
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/AdjointCondition/GetSecondDerivativesVector
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/ResponseFunction/CalculateGradient
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/ResponseFunction/CalculateFirstDerivativesGradient
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/strategies/strategies/explicit_strategy_testing_utils/AuxiliaryExplicitStrategiesTestElement/GetDofList
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/strategies/strategies/explicit_strategy_testing_utils/AuxiliaryExplicitStrategiesTestElement/AddExplicitContribution
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/strategies/strategies/explicit_strategy_testing_utils/AuxiliaryExplicitStrategiesTestElement/EquationIdVector
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/utilities/test_brute_force_point_locator/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/utilities/test_delaunator_utilities/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/utilities/test_discont_utils/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/utilities/test_divide_tetrahedra_3d_4/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/utilities/test_divide_triangle_2d_3/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/utilities/test_divide_triangle_3d_3/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/utilities/test_mls_shape_functions_utility/MLSShapeFunctionsUtilityConsistencyCheck
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/utilities/test_model_part_utils/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/utilities/test_sensitivity_builder/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_cad_json_input/TestCadJsonInput/test_cad_json_input_read
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_coordinate_transformation_utils/SetupModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_exact_integration/TestExactIntegration/test_line_exact_integration_1
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_exact_integration/TestExactIntegration/test_line_exact_integration_2
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_exact_integration/TestExactIntegration/test_line_exact_integration_3
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_exact_integration/TestExactIntegration/_test_triangle_exact_integration_1
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_exact_integration/TestExactIntegration/_test_triangle_exact_integration_2
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_exact_integration/TestExactIntegration/_test_triangle_exact_integration_3
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_exact_integration/TestExactIntegration/_test_quadrilateral_exact_integration_1
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_exact_integration/TestExactIntegration/_test_quadrilateral_exact_integration_2
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_geometry_data/TestGeometryData/test_geometry_data
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_model_part/TestModelPart/test_geometry_container
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_model_part_operation_utilities/TestModelPartOperationUtilities/test_Substract
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_model_part_utils_connectivity_generations/TestModelPartUtilsConnectivityGenerations/test_ModelPartUtils_FromConnectivityGenerateElements
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_model_part_utils_connectivity_generations/TestModelPartUtilsConnectivityGenerations/test_ModelPartUtils_FromConnectivityGenerateConditions
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_mortar_mapper/TestMortarMapperCore/__generate_discontinous_case
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_normal_utils/FiniteDifferenceNormalShapeSensitivityTest
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_obj_io/debug_vtk
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_obj_io/TestObjIO/test_WriteObjIO
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_processes/TestProcesses/test_assign_vector_variable_to_conditions_analytic_direction
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_processes/TestProcesses/test_split_internal_interfaces_process_3D
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_redistance/TestRedistance/test_parallel_redistance_maintain_plane_3d
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_sensitivity_utilities/CalculateEntityNeighbourMaps
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_sensitivity_utilities/TestSensitivityUtilitiesTwoDimSymmetricalSquare/__InitializeAssignEntityDerivativesToNodesValues
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_stl_io/TestStlIO/test_WriteStlIO
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_tensor_adaptors/TestTensorAdaptors/testGeometryMetricsTensorAdaptorDomainSizeCondition
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_tensor_adaptors/TestTensorAdaptors/testGeometryMetricsTensorAdaptorDomainSizeElement
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_utilities/test_bar_element/TestBarElement/Create
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_utilities/test_bar_element/TestBarElement/Clone
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_utilities/test_bar_element/TestBarElement/GetDofList
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_utilities/test_bar_element/TestBarElement/EquationIdVector
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_utilities/test_bar_element/TestBarElement/GetValuesVector
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_utilities/test_bar_element/TestBarElement/GetFirstDerivativesVector
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_utilities/test_bar_element/TestBarElement/GetSecondDerivativesVector
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_utilities/test_bar_element/TestBarElement/CalculateLeftHandSide
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_utilities/test_bar_element/TestBarElement/CalculateMassMatrix
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_utilities/test_bar_element/TestBarElement/CalculateLumpedMassVector
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_utilities/test_bar_element/TestBarElement/Check
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_utilities/test_element/TestElement/Create
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_utilities/test_element/TestElement/Clone
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_utilities/test_element/TestElement/GetDofList
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_utilities/test_element/TestElement/EquationIdVector
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_utilities/test_element/TestElement/GetValuesVector
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_utilities/test_element/TestElement/GetFirstDerivativesVector
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_utilities/test_element/TestElement/GetSecondDerivativesVector
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_utilities/test_element/TestElement/CalculateRightHandSide
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_utilities/test_element/TestElement/CalculateLeftHandSide
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_utilities/test_element/TestElement/CalculateMassMatrix
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_utilities/test_element/TestElement/CalculateDampingMatrix
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_utilities/test_element/TestElement/Check
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_utilities/test_element/TestElement/Initialize
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_utilities/test_laplacian_element/TestLaplacianElement/Create
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_utilities/test_laplacian_element/TestLaplacianElement/Clone
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_utilities/test_laplacian_element/TestLaplacianElement/GetDofList
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_utilities/test_laplacian_element/TestLaplacianElement/EquationIdVector
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_utilities/test_laplacian_element/TestLaplacianElement/GetValuesVector
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_utilities/test_laplacian_element/TestLaplacianElement/CalculateLocalSystem
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_utilities/test_laplacian_element/TestLaplacianElement/CalculateMassMatrix
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_utilities/test_laplacian_element/TestLaplacianElement/CalculateDampingMatrix
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_utilities/test_laplacian_element/TestLaplacianElement/Check
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_variable_redistribution/TestVariableRedistributionUtility/test_nodal_area
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_voxel_modeler/TestVoxelMeshModeler/test_TestVoxelMeshModelerNoMdpa
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_voxel_modeler/TestVoxelMeshModeler/test_TestVoxelMeshModelerElementInput
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_voxel_modeler/TestVoxelMeshModeler/test_TestVoxelMeshModelerQuadraticElementInput
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_voxel_modeler/TestVoxelMeshModeler/test_TestVoxelMeshModelerElementInputAutoSkin
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_voxel_modeler/TestVoxelMeshModeler/test_TestVoxelMeshModelerTetrahedralElementInput
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_voxel_modeler/TestVoxelMeshModeler/test_TestVoxelMeshModelerInRadius
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_voxel_modeler/TestBoundaryConditionVoxelizer/test_TestWallConditionVoxelizer
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_voxel_modeler/TestBoundaryConditionVoxelizer/test_TestFacesBetweenColors
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_voxel_modeler/TestBoundaryConditionVoxelizer/test_TestOuterWallConditionVoxelizer
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_vtu_output/fill_nodes
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_vtu_output/TestVtuOutput/GetNodes
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_vtu_output/TestVtuOutput/GetNumberOfNodes
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_vtu_output/TestVtuOutput/CheckGaussVtuFile
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/activation_utilities/ActivationUtilities/ActivateElementsAndConditions
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/auxiliar_model_part_utilities/AuxiliarModelPartUtilities/AddElementWithNodes
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/auxiliar_model_part_utilities/AuxiliarModelPartUtilities/AddConditionWithNodes
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/auxiliar_model_part_utilities/AuxiliarModelPartUtilities/RemoveElementAndBelongings
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/auxiliar_model_part_utilities/AuxiliarModelPartUtilities/RemoveElementsAndBelongings
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/auxiliar_model_part_utilities/AuxiliarModelPartUtilities/RemoveConditionAndBelongings
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/auxiliar_model_part_utilities/AuxiliarModelPartUtilities/RemoveConditionsAndBelongings
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/auxiliar_model_part_utilities/AuxiliarModelPartUtilities/RemoveOrphanNodesFromSubModelParts
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/auxiliar_model_part_utilities/AuxiliarModelPartUtilities/RetrieveElementsNeighbourElementsIds
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/auxiliar_model_part_utilities/AuxiliarModelPartUtilities/RetrieveConditionsNeighbourConditionsIds
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/auxiliar_model_part_utilities/AuxiliaryAddEntitiesWithNodes
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/binbased_fast_point_locator/FindPointOnMesh
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/binbased_nodes_in_element_locator/BinBasedNodesInElementLocator/FindNodesInElement
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/body_distance_calculation_utils/CalculateDistances
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/body_normal_calculation_utils/BodyNormalCalculationUtils/CalculateBodyNormals
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/brute_force_point_locator/BruteForcePointLocator/FindObject
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/convect_particles_utilities/ParticleConvectUtily/MoveParticles_Substepping
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/convect_particles_utilities/ParticleConvectUtily/MoveParticles_RK4
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/convect_particles_utilities/ParticleConvectUtily/EraseOuterElements
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/cutting_utility/CuttingUtility/FindSmallestEdge
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/cutting_utility/CuttingUtility/AddSkinConditions
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/cutting_utility/CuttingUtility/FirstLoop
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/cutting_utility/CuttingUtility/GenerateElements
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/embedded_skin_utility/EmbeddedSkinUtility<TDim>/GenerateSkin
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/embedded_skin_utility/EmbeddedSkinUtility<TDim>/ComputeElementSkin
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/embedded_skin_utility/EmbeddedSkinUtility<TDim>/SetDistancesVector
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/embedded_skin_utility/InterpolateMeshVariableToSkinSpecialization
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/entities_utilities/EntitityIdentifier<TEntity>/GenerateSingleType
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/entities_utilities/EntitityIdentifier<TEntity>/GenerateMultipleTypes
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/entities_utilities/EntitityIdentifier<TEntity>/GenerateTemplatedTypes
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/exact_mortar_segmentation_utility/ExactMortarIntegrationUtility<TDim, TNumNodes, TBelong, TNumNodesMaster>/TestGetExactIntegration
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/exact_mortar_segmentation_utility/ExactMortarIntegrationUtility<TDim, TNumNodes, TBelong, TNumNodesMaster>/TestGetExactAreaIntegration
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/exact_mortar_segmentation_utility/ExactMortarIntegrationUtility<TDim, TNumNodes, TBelong, TNumNodesMaster>/TestIODebug
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/input_output_utilities/InputOutputUtilities/SkippableEntity
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/iso_printer/IsosurfacePrinterApplication/AddSkinConditions
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/iso_printer/IsosurfacePrinterApplication/AddScalarVarIsosurface
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/iso_printer/IsosurfacePrinterApplication/AddScalarVarIsosurfaceAndLower
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/iso_printer/IsosurfacePrinterApplication/AddScalarVarIsosurfaceAndHigher
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/model_part_graph_utilities/ModelPartGraphUtilities/ComputeGraph
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/model_part_operation_utilities/ModelPartOperationUtilities/FillNodePointersForEntities
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/model_part_operation_utilities/ModelPartOperationUtilities/FillNodesFromEntities
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/model_part_operation_utilities/ModelPartOperationUtilities/FindNeighbourEntities
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/model_part_operation_utilities/ModelPartOperationUtilities/CheckEntities
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/model_part_operation_utilities/AddEntities
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/model_part_utils/GenerateEntitiesFromConnectivities
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/mortar_utilities/InvertNormalForFlag
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/normal_calculation_utils/NormalCalculationUtils/CalculateNormalsInContainer
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/normal_calculation_utils/NormalCalculationUtils/InitializeNormals
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/normal_calculation_utils/NormalCalculationUtils/CalculateNormalShapeDerivativesOnSimplex
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/normal_calculation_utils/for
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/normal_calculation_utils/NormalCalculationUtils/ComputeUnitNormalsFromAreaNormals
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/normal_calculation_utils/NormalCalculationUtils/CalculateNormal2D
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/normal_calculation_utils/NormalCalculationUtils/CalculateNormal3D
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/normal_calculation_utils/NormalCalculationUtils/CalculateNormalShapeDerivative3D
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/normal_calculation_utils/NormalCalculationUtils/CalculateNormalsUsingGenericAlgorithm
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/normal_calculation_utils/NormalCalculationUtils/CheckUseSimplex
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/normal_calculation_utils/NormalCalculationUtils/AuxiliaryCalculateOnSimplex
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/normal_calculation_utils/CalculateOnSimplex
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/normal_calculation_utils/CalculateOnSimplexLowMemory
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/particles_utilities/CountParticlesInNodes
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/particles_utilities/InterpolateValuesAtCoordinates
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/qef_utility/QuadraticErrorFunction/QuadraticErrorFunctionPoint
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/rve_periodicity_utility/RVEPeriodicityUtility/AuxiliaryAssignPeriodicity
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/sensitivity_builder/AssembleContainerContributions
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/sensitivity_utilities/SensitivityUtilities/AssignEntityDerivativesToNodes
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/shifted_boundary_meshless_interface_utility/ShiftedBoundaryMeshlessInterfaceUtility/CalculateGradientBasedConformingExtensionBasis
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/shifted_boundary_meshless_interface_utility/ShiftedBoundaryMeshlessInterfaceUtility/CalculateMeshlessBasedConformingExtensionBasis
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/shifted_boundary_meshless_interface_utility/ShiftedBoundaryMeshlessInterfaceUtility/CalculateMeshlessBasedNonConformingExtensionBasis
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/shifted_boundary_meshless_interface_utility/ShiftedBoundaryMeshlessInterfaceUtility/SetInterfaceFlags
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/shifted_boundary_meshless_interface_utility/ShiftedBoundaryMeshlessInterfaceUtility/SetSplitElementSupportCloud
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/shifted_boundary_meshless_interface_utility/ShiftedBoundaryMeshlessInterfaceUtility/SetSurrogateBoundaryNodalGradientWeights
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/signed_distance_calculation_utils/SignedDistanceCalculationUtils/CalculateDistances2D
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/signed_distance_calculation_utils/SignedDistanceCalculationUtils/CalculateDistances3D
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/spatial_containers_configure/SpatialContainersConfigure/CalculateBoundingBox
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/spatial_containers_configure/SpatialContainersConfigure/Intersection
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/spatial_containers_configure/SpatialContainersConfigure/IntersectionBox
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/specifications_utilities/SpecificationsUtilities/DetermineIfCompatibleGeometries
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/tessellation_utilities/delaunator_utilities/CreateTriangleMeshFromNodes
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/variable_redistribution_utility/VariableRedistributionUtility/CallSpecializedConvertDistributedValuesToPoint
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/variable_redistribution_utility/VariableRedistributionUtility/CallSpecializedDistributePointValues
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/variable_redistribution_utility/VariableRedistributionUtility/SpecializedConvertDistributedValuesToPoint
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/variable_redistribution_utility/VariableRedistributionUtility/ComputeNodalSizes
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/variable_redistribution_utility/VariableRedistributionUtility/UpdateDistributionRHS
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/variable_utils/VariableUtils/WeightedAccumulateVariableOnNodes
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static typename GeometryType::Pointer GetGeometry( const Parameters rParameters, ModelPart& rModelPart)`

# Calls

- [GetInt](../../../../functions/sources/kratos_parameters/Parameters/GetInt.md)
- [GetString](../../../../functions/sources/kratos_parameters/Parameters/GetString.md)

# Called by

- [Create](../../../../functions/conditions/mesh_condition/MeshCondition/Create.md)
- [Clone](../../../../functions/conditions/mesh_condition/MeshCondition/Clone.md)
- [CalculateOnIntegrationPoints](../../../../functions/conditions/mesh_condition/MeshCondition/CalculateOnIntegrationPoints.md)
- [CalculateEdgeDataValues](../../../../functions/containers/edge_based_data_structure/EdgeBasedDataStructure/CalculateEdgeDataValues.md)
- [Create](../../../../functions/elements/distance_calculation_element_simplex/DistanceCalculationElementSimplex/Create.md)
- [CalculateLocalSystem](../../../../functions/elements/distance_calculation_element_simplex/DistanceCalculationElementSimplex/CalculateLocalSystem.md)
- [EquationIdVector](../../../../functions/elements/distance_calculation_element_simplex/DistanceCalculationElementSimplex/EquationIdVector.md)
- [GetDofList](../../../../functions/elements/distance_calculation_element_simplex/DistanceCalculationElementSimplex/GetDofList.md)
- [Check](../../../../functions/elements/distance_calculation_element_simplex/DistanceCalculationElementSimplex/Check.md)
- [CalculateLocalSystem](../../../../functions/elements/edge_based_gradient_recovery_element/EdgeBasedGradientRecoveryElement<TDim>/CalculateLocalSystem.md)
- [Check](../../../../functions/elements/edge_based_gradient_recovery_element/EdgeBasedGradientRecoveryElement<TDim>/Check.md)
- [CalculateLeftHandSide](../../../../functions/elements/edge_based_gradient_recovery_element/EdgeBasedGradientRecoveryElement<TDim>/CalculateLeftHandSide.md)
- [CalculateRightHandSide](../../../../functions/elements/edge_based_gradient_recovery_element/EdgeBasedGradientRecoveryElement<TDim>/CalculateRightHandSide.md)
- [EquationIdVector](../../../../functions/elements/edge_based_gradient_recovery_element/EdgeBasedGradientRecoveryElement<2>/EquationIdVector.md)
- [EquationIdVector](../../../../functions/elements/edge_based_gradient_recovery_element/EdgeBasedGradientRecoveryElement<3>/EquationIdVector.md)
- [GetDofList](../../../../functions/elements/edge_based_gradient_recovery_element/EdgeBasedGradientRecoveryElement<2>/GetDofList.md)
- [GetDofList](../../../../functions/elements/edge_based_gradient_recovery_element/EdgeBasedGradientRecoveryElement<3>/GetDofList.md)
- [Create](../../../../functions/elements/edge_based_gradient_recovery_element/Create.md)
- [Check](../../../../functions/elements/embedded_nodal_variable_calculation_element_simplex/EmbeddedNodalVariableCalculationElementSimplex<double>/Check.md)
- [Check](../../../../functions/elements/embedded_nodal_variable_calculation_element_simplex/EmbeddedNodalVariableCalculationElementSimplex<array_1d<double,3>>/Check.md)
- [CalculateRightHandSide](../../../../functions/elements/embedded_nodal_variable_calculation_element_simplex/EmbeddedNodalVariableCalculationElementSimplex<double>/CalculateRightHandSide.md)
- [CalculateRightHandSide](../../../../functions/elements/embedded_nodal_variable_calculation_element_simplex/EmbeddedNodalVariableCalculationElementSimplex<array_1d<double, 3>>/CalculateRightHandSide.md)
- [EquationIdVector](../../../../functions/elements/embedded_nodal_variable_calculation_element_simplex/EmbeddedNodalVariableCalculationElementSimplex<double>/EquationIdVector.md)
- [EquationIdVector](../../../../functions/elements/embedded_nodal_variable_calculation_element_simplex/EmbeddedNodalVariableCalculationElementSimplex<array_1d<double, 3>>/EquationIdVector.md)
- [GetDofList](../../../../functions/elements/embedded_nodal_variable_calculation_element_simplex/EmbeddedNodalVariableCalculationElementSimplex<double>/GetDofList.md)
- [GetDofList](../../../../functions/elements/embedded_nodal_variable_calculation_element_simplex/EmbeddedNodalVariableCalculationElementSimplex<array_1d<double, 3>>/GetDofList.md)
- [Create](../../../../functions/elements/embedded_nodal_variable_calculation_element_simplex/Create.md)
- [Create](../../../../functions/elements/levelset_convection_element_simplex/LevelSetConvectionElementSimplex/Create.md)
- [Initialize](../../../../functions/elements/levelset_convection_element_simplex/LevelSetConvectionElementSimplex/Initialize.md)
- [CalculateLocalSystem](../../../../functions/elements/levelset_convection_element_simplex/LevelSetConvectionElementSimplex/CalculateLocalSystem.md)
- [EquationIdVector](../../../../functions/elements/levelset_convection_element_simplex/LevelSetConvectionElementSimplex/EquationIdVector.md)
- [GetDofList](../../../../functions/elements/levelset_convection_element_simplex/LevelSetConvectionElementSimplex/GetDofList.md)
- [Create](../../../../functions/elements/levelset_convection_element_simplex_algebraic_stabilization/LevelSetConvectionElementSimplexAlgebraicStabilization<TDim, TNumNodes>/Create.md)
- [CalculateLocalSystem](../../../../functions/elements/levelset_convection_element_simplex_algebraic_stabilization/LevelSetConvectionElementSimplexAlgebraicStabilization<TDim, TNumNodes>/CalculateLocalSystem.md)
- [Create](../../../../functions/elements/mesh_element/MeshElement/Create.md)
- [Clone](../../../../functions/elements/mesh_element/MeshElement/Clone.md)
- [CalculateOnIntegrationPoints](../../../../functions/elements/mesh_element/MeshElement/CalculateOnIntegrationPoints.md)
- [CalculateLeftHandSide](../../../../functions/future/tests/cpp_tests/test_utilities/solving_strategies_test_utilities/TestLaplacianElement1D/CalculateLeftHandSide.md)
- [CalculateRightHandSide](../../../../functions/future/tests/cpp_tests/test_utilities/solving_strategies_test_utilities/TestLaplacianElement1D/CalculateRightHandSide.md)
- [EquationIdVector](../../../../functions/future/tests/cpp_tests/test_utilities/solving_strategies_test_utilities/TestLaplacianElement1D/EquationIdVector.md)
- [GetDofList](../../../../functions/future/tests/cpp_tests/test_utilities/solving_strategies_test_utilities/TestLaplacianElement1D/GetDofList.md)
- [CalculateLeftHandSide](../../../../functions/future/tests/cpp_tests/test_utilities/solving_strategies_test_utilities/TestElasticityElement2D/CalculateLeftHandSide.md)
- [CalculateRightHandSide](../../../../functions/future/tests/cpp_tests/test_utilities/solving_strategies_test_utilities/TestElasticityElement2D/CalculateRightHandSide.md)
- [EquationIdVector](../../../../functions/future/tests/cpp_tests/test_utilities/solving_strategies_test_utilities/TestElasticityElement2D/EquationIdVector.md)
- [GetDofList](../../../../functions/future/tests/cpp_tests/test_utilities/solving_strategies_test_utilities/TestElasticityElement2D/GetDofList.md)
- [Create](../../../../functions/includes/condition/Condition/Create.md)
- [Clone](../../../../functions/includes/condition/Condition/Clone.md)
- [Check](../../../../functions/includes/condition/Condition/Check.md)
- [Create](../../../../functions/includes/element/Element/Create.md)
- [Clone](../../../../functions/includes/element/Element/Clone.md)
- [Check](../../../../functions/includes/element/Element/Check.md)
- [SetData](../../../../functions/includes/geometrical_object/SetData.md)
- [HasSameGeometryType](../../../../functions/includes/geometrical_object/HasSameGeometryType.md)
- [AddElement](../../../../functions/includes/gid_gauss_point_container/GidGaussPointsContainer/AddElement.md)
- [AddCondition](../../../../functions/includes/gid_gauss_point_container/GidGaussPointsContainer/AddCondition.md)
- [WriteSphereMesh](../../../../functions/includes/gid_io/WriteSphereMesh.md)
- [WriteClusterMesh](../../../../functions/includes/gid_io/WriteClusterMesh.md)
- [AddElement](../../../../functions/includes/gid_mesh_container/GidMeshContainer/AddElement.md)
- [AddCondition](../../../../functions/includes/gid_mesh_container/GidMeshContainer/AddCondition.md)
- [WriteMesh](../../../../functions/includes/gid_mesh_container/GidMeshContainer/WriteMesh.md)
- [WorkingSpaceDimension](../../../../functions/includes/mesh/Mesh/WorkingSpaceDimension.md)
- [ReadBrepEdgeBrepCurveOnSurface](../../../../functions/input_output/cad_json_input/CadJsonInput/ReadBrepEdgeBrepCurveOnSurface.md)
- [ReadCouplingGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/ReadCouplingGeometry.md)
- [ReadPointsOnGeometries](../../../../functions/input_output/cad_json_input/CadJsonInput/ReadPointsOnGeometries.md)
- [WriteGeometricalGaussVariableToFile](../../../../functions/input_output/ensight_output/EnSightOutput/WriteGeometricalGaussVariableToFile.md)
- [CollectPartData](../../../../functions/input_output/ensight_output/EnSightOutput/CollectPartData.md)
- [if](../../../../functions/input_output/ensight_output/if.md)
- [GetGeometryConnectivity](../../../../functions/input_output/ensight_output/EnSightOutput/GetGeometryConnectivity.md)
- [GetAverageIntegrationValue](../../../../functions/input_output/ensight_output/EnSightOutput/GetAverageIntegrationValue.md)
- [WriteModelPart](../../../../functions/input_output/obj_io/ObjIO/WriteModelPart.md)
- [WriteEntityBlock](../../../../functions/input_output/stl_io/StlIO/WriteEntityBlock.md)
- [WriteEntityBlockMPI](../../../../functions/input_output/stl_io/StlIO/WriteEntityBlockMPI.md)
- [WriteElements](../../../../functions/input_output/unv_output/UnvOutput/WriteElements.md)
- [DetermineVtkCellListSize](../../../../functions/input_output/vtk_output/VtkOutput/DetermineVtkCellListSize.md)
- [WriteCellType](../../../../functions/input_output/vtk_output/VtkOutput/WriteCellType.md)
- [WriteIntegrationScalarContainerVariable](../../../../functions/input_output/vtk_output/VtkOutput/WriteIntegrationScalarContainerVariable.md)
- [WriteIntegrationVectorContainerVariable](../../../../functions/input_output/vtk_output/VtkOutput/WriteIntegrationVectorContainerVariable.md)
- [WriteModelPartWithoutNodesToFile](../../../../functions/input_output/vtk_output/VtkOutput/WriteModelPartWithoutNodesToFile.md)
- [GetGeometryTypes](../../../../functions/input_output/vtu_output/GetGeometryTypes.md)
- [GetOffsets](../../../../functions/input_output/vtu_output/GetOffsets.md)
- [GetConnectivities](../../../../functions/input_output/vtu_output/GetConnectivities.md)
- [GetNodesContainer](../../../../functions/input_output/vtu_output/GetNodesContainer.md)
- [WriteIntegrationPointData](../../../../functions/input_output/vtu_output/VtuOutput/WriteIntegrationPointData.md)
- [GenerateMesh](../../../../functions/modeler/cartesian_mesh_generator_modeler/CartesianMeshGeneratorModeler/GenerateMesh.md)
- [CalculateNormals](../../../../functions/modeler/cartesian_mesh_generator_modeler/CartesianMeshGeneratorModeler/CalculateNormals.md)
- [CalculateBoundaryIntersections](../../../../functions/modeler/cartesian_mesh_generator_modeler/CartesianMeshGeneratorModeler/CalculateBoundaryIntersections.md)
- [CalculateBoundingBox](../../../../functions/modeler/cartesian_mesh_generator_modeler/CartesianMeshGeneratorModeler/CalculateBoundingBox.md)
- [CleanUpProblematicGeometries](../../../../functions/modeler/clean_up_problematic_triangles_modeler/CleanUpProblematicTrianglesModeler/CleanUpProblematicGeometries.md)
- [ComputeNullAreaTriangles](../../../../functions/modeler/clean_up_problematic_triangles_modeler/CleanUpProblematicTrianglesModeler/ComputeNullAreaTriangles.md)
- [Apply](../../../../functions/modeler/coloring/color_and_intersect_with_cells_in_touch/ColorAndIntersectWithCellsInTouch/Apply.md)
- [Apply](../../../../functions/modeler/coloring/color_cell_faces/ColorCellFaces/Apply.md)
- [Apply](../../../../functions/modeler/coloring/color_cells_in_touch/ColorCellsInTouch/Apply.md)
- [Apply](../../../../functions/modeler/coloring/color_cells_with_inside_center/ColorCellsWithInsideCenter/Apply.md)
- [Apply](../../../../functions/modeler/coloring/color_connected_cells_in_touch/ColorConnectedCellsInTouch/Apply.md)
- [GenerateMesh](../../../../functions/modeler/duplicate_mesh_modeler/DuplicateMeshModeler/GenerateMesh.md)
- [WriteMesh](../../../../functions/modeler/edge_swapping_2d_modeler/EdgeSwapping2DModeler/WriteMesh.md)
- [CollapseNodes](../../../../functions/modeler/edge_swapping_2d_modeler/EdgeSwapping2DModeler/CollapseNodes.md)
- [SetNodalCollapsingData](../../../../functions/modeler/edge_swapping_2d_modeler/EdgeSwapping2DModeler/SetNodalCollapsingData.md)
- [Swap](../../../../functions/modeler/edge_swapping_2d_modeler/EdgeSwapping2DModeler/Swap.md)
- [ElementQuality](../../../../functions/modeler/edge_swapping_2d_modeler/EdgeSwapping2DModeler/ElementQuality.md)
- [FindNodalNeighbours](../../../../functions/modeler/edge_swapping_2d_modeler/EdgeSwapping2DModeler/FindNodalNeighbours.md)
- [SetSwappingData](../../../../functions/modeler/edge_swapping_2d_modeler/EdgeSwapping2DModeler/SetSwappingData.md)
- [FindNeighbourElement](../../../../functions/modeler/edge_swapping_2d_modeler/EdgeSwapping2DModeler/FindNeighbourElement.md)
- [IsInElementSphere](../../../../functions/modeler/edge_swapping_2d_modeler/EdgeSwapping2DModeler/IsInElementSphere.md)
- [GenerateIntersection](../../../../functions/modeler/entity_generation/generate_stl_intersection_with_cells/GenerateStlIntersectionWithCells/GenerateIntersection.md)
- [SetMesh](../../../../functions/modeler/mesh_suite_modeler/MeshSuiteModeler/SetMesh.md)
- [Execute](../../../../functions/modeler/operation/compute_surrogate_boundary_data/ComputeSurrogateBoundaryData/Execute.md)
- [FindConditionContact](../../../../functions/modeler/operation/find_contacts_in_skin_model_part/FindContactsInSkinModelPart/FindConditionContact.md)
- [GenerateMesh](../../../../functions/modeler/structured_mesh_refinement_modeler/StructuredMeshRefinementModeler/GenerateMesh.md)
- [GenerateNodes](../../../../functions/modeler/structured_mesh_refinement_modeler/GenerateNodes.md)
- [GenerateNode](../../../../functions/modeler/structured_mesh_refinement_modeler/GenerateNode.md)
- [GenerateElements](../../../../functions/modeler/structured_mesh_refinement_modeler/StructuredMeshRefinementModeler/GenerateElements.md)
- [GenerateConditions](../../../../functions/modeler/structured_mesh_refinement_modeler/StructuredMeshRefinementModeler/GenerateConditions.md)
- [__InitializeModelPart](../../../../functions/mpi/future/tests/test_mpi_vtu_output/TestMPIVtuOutput/__InitializeModelPart.md)
- [executeComputeArea_Task](../../../../functions/mpi/tests/test_mpi_serializer/executeComputeArea_Task.md)
- [GatherModelPartUtility](../../../../functions/mpi/utilities/gather_modelpart_utility/GatherModelPartUtility/GatherModelPartUtility.md)
- [OrientFaces](../../../../functions/mpi/utilities/mpi_normal_calculation_utilities/MPINormalCalculationUtils/OrientFaces.md)
- [IdentifyFaces](../../../../functions/mpi/utilities/mpi_normal_calculation_utilities/MPINormalCalculationUtils/IdentifyFaces.md)
- [InitializeNormalData](../../../../functions/mpi/utilities/mpi_normal_calculation_utilities/MPINormalCalculationUtils/InitializeNormalData.md)
- [UpdateNodeNormals](../../../../functions/mpi/utilities/mpi_normal_calculation_utilities/MPINormalCalculationUtils/UpdateNodeNormals.md)
- [ApplyConstraintsForPeriodicConditions](../../../../functions/processes/apply_periodic_boundary_condition_process/ApplyPeriodicConditionProcess/ApplyConstraintsForPeriodicConditions.md)
- [GetCellIntersections](../../../../functions/processes/apply_ray_casting_process/ApplyRayCastingProcess<TDim>/GetCellIntersections.md)
- [CallFunction](../../../../functions/processes/assign_scalar_field_to_entities_process/AssignScalarFieldToEntitiesProcess<Condition>/CallFunction.md)
- [CallFunction](../../../../functions/processes/assign_scalar_field_to_entities_process/AssignScalarFieldToEntitiesProcess<Element>/CallFunction.md)
- [CallFunctionComponents](../../../../functions/processes/assign_scalar_field_to_entities_process/AssignScalarFieldToEntitiesProcess<Condition>/CallFunctionComponents.md)
- [CallFunctionComponents](../../../../functions/processes/assign_scalar_field_to_entities_process/AssignScalarFieldToEntitiesProcess<Element>/CallFunctionComponents.md)
- [CallFunctionLocalSystem](../../../../functions/processes/assign_scalar_field_to_entities_process/AssignScalarFieldToEntitiesProcess<Condition>/CallFunctionLocalSystem.md)
- [CallFunctionLocalSystem](../../../../functions/processes/assign_scalar_field_to_entities_process/AssignScalarFieldToEntitiesProcess<Element>/CallFunctionLocalSystem.md)
- [CallFunctionLocalSystemComponents](../../../../functions/processes/assign_scalar_field_to_entities_process/AssignScalarFieldToEntitiesProcess<Condition>/CallFunctionLocalSystemComponents.md)
- [CallFunctionLocalSystemComponents](../../../../functions/processes/assign_scalar_field_to_entities_process/AssignScalarFieldToEntitiesProcess<Element>/CallFunctionLocalSystemComponents.md)
- [AssignTimeDependentValue](../../../../functions/processes/assign_scalar_field_to_entities_process/AssignScalarFieldToEntitiesProcess<Condition>/AssignTimeDependentValue.md)
- [AssignTimeDependentValue](../../../../functions/processes/assign_scalar_field_to_entities_process/AssignScalarFieldToEntitiesProcess<Element>/AssignTimeDependentValue.md)
- [GetCoordinatesEntity](../../../../functions/processes/assign_scalar_input_to_entities_process/AssignScalarInputToEntitiesProcess<Condition, AssignScalarInputToEntitiesProcessSettings/SaveAsNonHistoricalVariable>/GetCoordinatesEntity.md)
- [GetCoordinatesEntity](../../../../functions/processes/assign_scalar_input_to_entities_process/AssignScalarInputToEntitiesProcess<Element, AssignScalarInputToEntitiesProcessSettings/SaveAsNonHistoricalVariable>/GetCoordinatesEntity.md)
- [CalculateElementalDistances](../../../../functions/processes/calculate_discontinuous_distance_to_skin_process/CalculateDiscontinuousDistanceToSkinProcess<TDim>/CalculateElementalDistances.md)
- [CalculateElementalAndEdgeDistances](../../../../functions/processes/calculate_discontinuous_distance_to_skin_process/CalculateDiscontinuousDistanceToSkinProcess<TDim>/CalculateElementalAndEdgeDistances.md)
- [ComputeEdgesIntersections](../../../../functions/processes/calculate_discontinuous_distance_to_skin_process/CalculateDiscontinuousDistanceToSkinProcess<TDim>/ComputeEdgesIntersections.md)
- [ComputeElementalDistancesFromPlaneApproximation](../../../../functions/processes/calculate_discontinuous_distance_to_skin_process/CalculateDiscontinuousDistanceToSkinProcess<TDim>/ComputeElementalDistancesFromPlaneApproximation.md)
- [ComputeIntersectionPlaneElementalDistances](../../../../functions/processes/calculate_discontinuous_distance_to_skin_process/CalculateDiscontinuousDistanceToSkinProcess<TDim>/ComputeIntersectionPlaneElementalDistances.md)
- [ComputePlaneApproximation](../../../../functions/processes/calculate_discontinuous_distance_to_skin_process/CalculateDiscontinuousDistanceToSkinProcess<TDim>/ComputePlaneApproximation.md)
- [CorrectDistanceOrientation](../../../../functions/processes/calculate_discontinuous_distance_to_skin_process/CalculateDiscontinuousDistanceToSkinProcess<TDim>/CorrectDistanceOrientation.md)
- [ComputeElementalDistancesFromEdgeRatios](../../../../functions/processes/calculate_discontinuous_distance_to_skin_process/CalculateDiscontinuousDistanceToSkinProcess<TDim>/ComputeElementalDistancesFromEdgeRatios.md)
- [CheckAndCorrectEdgeDistances](../../../../functions/processes/calculate_discontinuous_distance_to_skin_process/CalculateDiscontinuousDistanceToSkinProcess<TDim>/CheckAndCorrectEdgeDistances.md)
- [CalculateEmbeddedVariableFromSkinSpecialization](../../../../functions/processes/calculate_discontinuous_distance_to_skin_process/CalculateEmbeddedVariableFromSkinSpecialization.md)
- [Execute](../../../../functions/processes/calculate_distance_to_path_process/CalculateDistanceToPathProcess<THistorical>/Execute.md)
- [CalculateElementalDistances](../../../../functions/processes/calculate_distance_to_skin_process/CalculateDistanceToSkinProcess<TDim>/CalculateElementalDistances.md)
- [CalculateDistanceToNode](../../../../functions/processes/calculate_distance_to_skin_process/CalculateDistanceToSkinProcess<TDim>/CalculateDistanceToNode.md)
- [CalculateNodalDistances](../../../../functions/processes/calculate_distance_to_skin_process/CalculateDistanceToSkinProcess<TDim>/CalculateNodalDistances.md)
- [CalculateEmbeddedNodalVariableFromSkinProcess](../../../../functions/processes/calculate_embedded_nodal_variable_from_skin_process/CalculateEmbeddedNodalVariableFromSkinProcess/CalculateEmbeddedNodalVariableFromSkinProcess.md)
- [AddIntersectedElementsModelPartElements](../../../../functions/processes/calculate_embedded_nodal_variable_from_skin_process/CalculateEmbeddedNodalVariableFromSkinProcess/AddIntersectedElementsModelPartElements.md)
- [CalculateEmbeddedVelocity](../../../../functions/processes/calculate_embedded_signed_distance_to_3d_skin_process/CalculateEmbeddedSignedDistanceTo3DSkinProcess/CalculateEmbeddedVelocity.md)
- [SetMaximumAndMinimumDistanceValues](../../../../functions/processes/calculate_embedded_signed_distance_to_3d_skin_process/CalculateEmbeddedSignedDistanceTo3DSkinProcess/SetMaximumAndMinimumDistanceValues.md)
- [Execute](../../../../functions/processes/calculate_nodal_area_process/CalculateNodalAreaProcess<THistorical>/Execute.md)
- [CalculateBoundingBox](../../../../functions/processes/calculate_signed_distance_to_3d_condition_skin_process/DistanceSpatialContainersConditionConfigure/CalculateBoundingBox.md)
- [GetBoundingBox](../../../../functions/processes/calculate_signed_distance_to_3d_condition_skin_process/DistanceSpatialContainersConditionConfigure/GetBoundingBox.md)
- [Intersection](../../../../functions/processes/calculate_signed_distance_to_3d_condition_skin_process/DistanceSpatialContainersConditionConfigure/Intersection.md)
- [IntersectionBox](../../../../functions/processes/calculate_signed_distance_to_3d_condition_skin_process/DistanceSpatialContainersConditionConfigure/IntersectionBox.md)
- [IdentifyIntersectionNodes](../../../../functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/IdentifyIntersectionNodes.md)
- [CalcSignedDistancesToOneIntNode](../../../../functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/CalcSignedDistancesToOneIntNode.md)
- [CalcSignedDistancesToTwoIntNodes](../../../../functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/CalcSignedDistancesToTwoIntNodes.md)
- [CalcSignedDistancesToApproxTriangle](../../../../functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/CalcSignedDistancesToApproxTriangle.md)
- [AssignDistancesToElements](../../../../functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/AssignDistancesToElements.md)
- [CalculateNotEmptyLeavesDistance](../../../../functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/CalculateNotEmptyLeavesDistance.md)
- [GetCellIntersections](../../../../functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/GetCellIntersections.md)
- [CalculateBoundingBox](../../../../functions/processes/calculate_signed_distance_to_3d_skin_process/DistanceSpatialContainersConfigure/CalculateBoundingBox.md)
- [GetBoundingBox](../../../../functions/processes/calculate_signed_distance_to_3d_skin_process/DistanceSpatialContainersConfigure/GetBoundingBox.md)
- [Intersection](../../../../functions/processes/calculate_signed_distance_to_3d_skin_process/DistanceSpatialContainersConfigure/Intersection.md)
- [IntersectionBox](../../../../functions/processes/calculate_signed_distance_to_3d_skin_process/DistanceSpatialContainersConfigure/IntersectionBox.md)
- [MappingPressureToStructure](../../../../functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/MappingPressureToStructure.md)
- [AveragePressureToNode](../../../../functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/AveragePressureToNode.md)
- [IdentifyIntersectionNodes](../../../../functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/IdentifyIntersectionNodes.md)
- [CalcSignedDistancesToOneIntNode](../../../../functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/CalcSignedDistancesToOneIntNode.md)
- [CalcSignedDistancesToTwoIntNodes](../../../../functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/CalcSignedDistancesToTwoIntNodes.md)
- [CalcSignedDistancesToThreeIntNodes](../../../../functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/CalcSignedDistancesToThreeIntNodes.md)
- [CalcSignedDistancesToMoreThanThreeIntNodes](../../../../functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/CalcSignedDistancesToMoreThanThreeIntNodes.md)
- [AssignMinimalNodalDistance](../../../../functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/AssignMinimalNodalDistance.md)
- [GenerateSkinModelPart](../../../../functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/GenerateSkinModelPart.md)
- [CalculateNotEmptyLeavesDistance](../../../../functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/CalculateNotEmptyLeavesDistance.md)
- [GetCellIntersections](../../../../functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/GetCellIntersections.md)
- [Execute](../../../../functions/processes/check_same_modelpart_using_skin_distance_process/CheckSameModelPartUsingSkinDistanceProcess<TDim>/Execute.md)
- [Execute](../../../../functions/processes/check_skin_process/CheckSkinProcess/Execute.md)
- [ComputeElementalContributionsAndVolume](../../../../functions/processes/compute_nodal_gradient_process/ComputeNodalGradientProcess<TOutputHistorical>/ComputeElementalContributionsAndVolume.md)
- [Execute](../../../../functions/processes/compute_nodal_normal_divergence_process/ComputeNodalNormalDivergenceProcess<THistorical>/Execute.md)
- [CheckAndAssignSettings](../../../../functions/processes/edge_based_gradient_recovery_process/EdgeBasedGradientRecoveryProcess/CheckAndAssignSettings.md)
- [ReplicateWithoutFlags](../../../../functions/processes/fast_transfer_between_model_parts_process/FastTransferBetweenModelPartsProcess/ReplicateWithoutFlags.md)
- [ReplicateWithFlags](../../../../functions/processes/fast_transfer_between_model_parts_process/FastTransferBetweenModelPartsProcess/ReplicateWithFlags.md)
- [Execute](../../../../functions/processes/find_conditions_neighbours_process/FindConditionsNeighboursProcess/Execute.md)
- [CheckForNeighbourFaces](../../../../functions/processes/find_conditions_neighbours_process/FindConditionsNeighboursProcess/CheckForNeighbourFaces.md)
- [ComputeDimension](../../../../functions/processes/find_conditions_neighbours_process/FindConditionsNeighboursProcess/ComputeDimension.md)
- [Execute](../../../../functions/processes/find_global_nodal_entity_neighbours_process/FindGlobalNodalEntityNeighboursProcess<TContainerType>/Execute.md)
- [Execute](../../../../functions/processes/find_global_nodal_neighbours_for_entities_process/FindNodalNeighboursForEntitiesProcess<TContainerType>/Execute.md)
- [WorkingSpaceDimension](../../../../functions/processes/find_intersected_geometrical_objects_process/FindIntersectedGeometricalObjectsProcess/WorkingSpaceDimension.md)
- [MarkIfIntersected](../../../../functions/processes/find_intersected_geometrical_objects_process/FindIntersectedGeometricalObjectsProcess/MarkIfIntersected.md)
- [FindIntersectedSkinObjects](../../../../functions/processes/find_intersected_geometrical_objects_process/FindIntersectedGeometricalObjectsProcess/FindIntersectedSkinObjects.md)
- [CalculateBoundingBox](../../../../functions/processes/find_intersected_geometrical_objects_process/DistanceSpatialContainersConfigure/CalculateBoundingBox.md)
- [GetBoundingBox](../../../../functions/processes/find_intersected_geometrical_objects_process/DistanceSpatialContainersConfigure/GetBoundingBox.md)
- [Intersection](../../../../functions/processes/find_intersected_geometrical_objects_process/DistanceSpatialContainersConfigure/Intersection.md)
- [IntersectionBox](../../../../functions/processes/find_intersected_geometrical_objects_process/DistanceSpatialContainersConfigure/IntersectionBox.md)
- [Execute](../../../../functions/processes/find_nodal_h_process/FindNodalHProcess<THistorical>/Execute.md)
- [Execute](../../../../functions/processes/generic_find_elements_neighbours_process/GenericFindElementalNeighboursProcess/Execute.md)
- [ExecuteFinalizeSolutionStep](../../../../functions/processes/integration_values_extrapolation_to_nodes_process/IntegrationValuesExtrapolationToNodesProcess/ExecuteFinalizeSolutionStep.md)
- [InitializeMaps](../../../../functions/processes/integration_values_extrapolation_to_nodes_process/IntegrationValuesExtrapolationToNodesProcess/InitializeMaps.md)
- [EvaluateNumberOfSubsteps](../../../../functions/processes/levelset_convection_process/LevelSetConvectionProcess/EvaluateNumberOfSubsteps.md)
- [EvaluateLimiter](../../../../functions/processes/levelset_convection_process/LevelSetConvectionProcess/EvaluateLimiter.md)
- [InitializeConvectionStrategy](../../../../functions/processes/levelset_convection_process/LevelSetConvectionProcess/InitializeConvectionStrategy.md)
- [Execute](../../../../functions/processes/measure_mesh_quality_process/MeasureMeshQualityProcess/Execute.md)
- [SelectLowQualityElementNodes](../../../../functions/processes/mesh_local_smoothing_process/MeshLocalSmoothingProcess/SelectLowQualityElementNodes.md)
- [CollapseNode](../../../../functions/processes/mesh_node_collapsing_process/MeshNodeCollapsingProcess/CollapseNode.md)
- [CalculateMinQualityOfNeighbourElements](../../../../functions/processes/mesh_node_collapsing_process/MeshNodeCollapsingProcess/CalculateMinQualityOfNeighbourElements.md)
- [ElementHas](../../../../functions/processes/mesh_node_collapsing_process/MeshNodeCollapsingProcess/ElementHas.md)
- [SwapElementNode](../../../../functions/processes/mesh_node_collapsing_process/MeshNodeCollapsingProcess/SwapElementNode.md)
- [FindMaximumEdgeSize](../../../../functions/processes/parallel_distance_calculation_process/ParallelDistanceCalculationProcess<TDim>/FindMaximumEdgeSize.md)
- [SetDistancesOnDividedElements](../../../../functions/processes/parallel_distance_calculation_process/ParallelDistanceCalculationProcess<TDim>/SetDistancesOnDividedElements.md)
- [ExtendDistancesByLayer](../../../../functions/processes/parallel_distance_calculation_process/ParallelDistanceCalculationProcess<TDim>/ExtendDistancesByLayer.md)
- [Execute](../../../../functions/processes/reorder_and_optimize_modelpart_process/ReorderAndOptimizeModelPartProcess/Execute.md)
- [OptimizeOrdering](../../../../functions/processes/reorder_and_optimize_modelpart_process/ReorderAndOptimizeModelPartProcess/OptimizeOrdering.md)
- [ReorderElements](../../../../functions/processes/reorder_and_optimize_modelpart_process/ReorderAndOptimizeModelPartProcess/ReorderElements.md)
- [ExecuteInitializeSolutionStep](../../../../functions/processes/set_initial_state_process/SetInitialStateProcess<TDim>/ExecuteInitializeSolutionStep.md)
- [GetReferenceArea](../../../../functions/processes/simple_mortar_mapper_process/SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>/GetReferenceArea.md)
- [PerformMortarOperations](../../../../functions/processes/simple_mortar_mapper_process/PerformMortarOperations.md)
- [ClearIndexes](../../../../functions/processes/simple_mortar_mapper_process/ClearIndexes.md)
- [FillDatabase](../../../../functions/processes/simple_mortar_mapper_process/FillDatabase.md)
- [SimpleMortarMapperProcessWrapper](../../../../functions/processes/simple_mortar_mapper_wrapper_process/SimpleMortarMapperProcessWrapper/SimpleMortarMapperProcessWrapper.md)
- [GenerateFaceMaps](../../../../functions/processes/skin_detection_process/SkinDetectionProcess<TDim>/GenerateFaceMaps.md)
- [SetUpAdditionalSubModelParts](../../../../functions/processes/skin_detection_process/SkinDetectionProcess<TDim>/SetUpAdditionalSubModelParts.md)
- [ExecuteInitialize](../../../../functions/processes/split_internal_interfaces_process/SplitInternalInterfacesProcess/ExecuteInitialize.md)
- [SplitBoundary](../../../../functions/processes/split_internal_interfaces_process/SplitInternalInterfacesProcess/SplitBoundary.md)
- [Execute](../../../../functions/processes/tetrahedra_mesh_edge_swapping_process/TetrahedraMeshEdgeSwappingProcess/Execute.md)
- [FindOptimumPositionsAndWeights](../../../../functions/processes/tetrahedra_mesh_quality_weighted_smoothing_process/TetrahedraMeshQualityWeightedSmoothingProcess/FindOptimumPositionsAndWeights.md)
- [FindOptimumPositionsAndWeights](../../../../functions/processes/tetrahedra_mesh_worst_element_smoothing_process/TetrahedraMeshWorstElementSmoothingProcess/FindOptimumPositionsAndWeights.md)
- [Execute](../../../../functions/processes/tetrahedral_mesh_orientation_check/TetrahedralMeshOrientationCheck/Execute.md)
- [SwapAll](../../../../functions/processes/tetrahedral_mesh_orientation_check/TetrahedralMeshOrientationCheck/SwapAll.md)
- [SwapNegativeElements](../../../../functions/processes/tetrahedral_mesh_orientation_check/TetrahedralMeshOrientationCheck/SwapNegativeElements.md)
- [Execute](../../../../functions/processes/variational_distance_calculation_process/VariationalDistanceCalculationProcess/Execute.md)
- [ValidateInput](../../../../functions/processes/variational_distance_calculation_process/VariationalDistanceCalculationProcess/ValidateInput.md)
- [ReGenerateDistanceModelPart](../../../../functions/processes/variational_distance_calculation_process/VariationalDistanceCalculationProcess/ReGenerateDistanceModelPart.md)
- [VectorizedFindHelper](../../../../functions/python/add_geometrical_utilities_to_python/VectorizedFindHelper.md)
- [GetNodeFromElement](../../../../functions/python/add_mesh_to_python/GetNodeFromElement.md)
- [GetNodesFromElement](../../../../functions/python/add_mesh_to_python/GetNodesFromElement.md)
- [GetNodeFromCondition](../../../../functions/python/add_mesh_to_python/GetNodeFromCondition.md)
- [GetNodesFromCondition](../../../../functions/python/add_mesh_to_python/GetNodesFromCondition.md)
- [GetIntegrationPointsFromElement](../../../../functions/python/add_mesh_to_python/GetIntegrationPointsFromElement.md)
- [SetValuesOnIntegrationPoints](../../../../functions/python/add_mesh_to_python/SetValuesOnIntegrationPoints.md)
- [SetValuesOnIntegrationPointsVector](../../../../functions/python/add_mesh_to_python/SetValuesOnIntegrationPointsVector.md)
- [SetValuesOnIntegrationPointsConstitutiveLaw](../../../../functions/python/add_mesh_to_python/SetValuesOnIntegrationPointsConstitutiveLaw.md)
- [MakePRestrictionOperator](../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/p_multigrid_utilities/MakePRestrictionOperator.md)
- [Build](../../../../functions/solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_componentwise/ResidualBasedEliminationBuilderAndSolverComponentwise/Build.md)
- [CalculatePreviousTimeStepContributions](../../../../functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/ResidualBasedAdjointBossakScheme/CalculatePreviousTimeStepContributions.md)
- [CalculateNodeNeighbourCount](../../../../functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/ResidualBasedAdjointBossakScheme/CalculateNodeNeighbourCount.md)
- [UpdateEntityTimeSchemeContributions](../../../../functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/UpdateEntityTimeSchemeContributions.md)
- [UpdateEntityAuxiliaryVariableContributions](../../../../functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/UpdateEntityAuxiliaryVariableContributions.md)
- [CalculateSystemContributions](../../../../functions/solving_strategies/schemes/residualbased_incrementalupdate_static_scheme_slip/ResidualBasedIncrementalUpdateStaticSchemeSlip/CalculateSystemContributions.md)
- [CalculateRHSContribution](../../../../functions/solving_strategies/schemes/residualbased_incrementalupdate_static_scheme_slip/ResidualBasedIncrementalUpdateStaticSchemeSlip/CalculateRHSContribution.md)
- [CalculateLHSContribution](../../../../functions/solving_strategies/schemes/residualbased_incrementalupdate_static_scheme_slip/ResidualBasedIncrementalUpdateStaticSchemeSlip/CalculateLHSContribution.md)
- [CalculateLocalSensitivityAndGlobalPointersVector](../../../../functions/solving_strategies/schemes/sensitivity_builder_scheme/CalculateLocalSensitivityAndGlobalPointersVector.md)
- [WriteElements](../../../../functions/sources/model_part_io/ModelPartIO/WriteElements.md)
- [WriteConditions](../../../../functions/sources/model_part_io/ModelPartIO/WriteConditions.md)
- [FillNodalConnectivitiesFromElementBlockInList](../../../../functions/sources/model_part_io/ModelPartIO/FillNodalConnectivitiesFromElementBlockInList.md)
- [FillNodalConnectivitiesFromConditionBlockInList](../../../../functions/sources/model_part_io/ModelPartIO/FillNodalConnectivitiesFromConditionBlockInList.md)
- [ReadElementsBlock](../../../../functions/sources/model_part_io/ModelPartIO/ReadElementsBlock.md)
- [ReadConditionsBlock](../../../../functions/sources/model_part_io/ModelPartIO/ReadConditionsBlock.md)
- [ReadElementsConnectivitiesBlock](../../../../functions/sources/model_part_io/ModelPartIO/ReadElementsConnectivitiesBlock.md)
- [ReadConditionsConnectivitiesBlock](../../../../functions/sources/model_part_io/ModelPartIO/ReadConditionsConnectivitiesBlock.md)
- [FillNodalConnectivitiesFromElementBlock](../../../../functions/sources/model_part_io/ModelPartIO/FillNodalConnectivitiesFromElementBlock.md)
- [FillNodalConnectivitiesFromConditionBlock](../../../../functions/sources/model_part_io/ModelPartIO/FillNodalConnectivitiesFromConditionBlock.md)
- [DivideElementsBlock](../../../../functions/sources/model_part_io/ModelPartIO/DivideElementsBlock.md)
- [DivideConditionsBlock](../../../../functions/sources/model_part_io/ModelPartIO/DivideConditionsBlock.md)
- [Create](../../../../functions/sources/periodic_condition/PeriodicCondition/Create.md)
- [GetDofList](../../../../functions/sources/periodic_condition/PeriodicCondition/GetDofList.md)
- [Partition](../../../../functions/spatial_containers/bounding_volume_tree/LineRegressionVolumePartitioner/Partition.md)
- [Partition](../../../../functions/spatial_containers/bounding_volume_tree/SimpleBoundingVolumePartitioner/Partition.md)
- [UpdateTree](../../../../functions/spatial_containers/bounding_volume_tree/BoundingVolumeTree/UpdateTree.md)
- [AddGeometryFromConditions](../../../../functions/spatial_containers/bounding_volume_tree/BoundingVolumeTree/AddGeometryFromConditions.md)
- [Initialize](../../../../functions/spatial_containers/dem_search/RadiusPoint/Initialize.md)
- [SearchInRadius](../../../../functions/spatial_containers/geometrical_objects_bins/GeometricalObjectsBins/SearchInRadius.md)
- [SearchNearestInCell](../../../../functions/spatial_containers/geometrical_objects_bins/GeometricalObjectsBins/SearchNearestInCell.md)
- [SearchIsInsideInCell](../../../../functions/spatial_containers/geometrical_objects_bins/GeometricalObjectsBins/SearchIsInsideInCell.md)
- [GeometricalObjectsBins](../../../../functions/spatial_containers/geometrical_objects_bins/GeometricalObjectsBins.md)
- [AddObjectsToCells](../../../../functions/spatial_containers/geometrical_objects_bins/AddObjectsToCells.md)
- [IsIntersected](../../../../functions/spatial_containers/octree_binary/OctreeBinary/IsIntersected.md)
- [UpdatePoint](../../../../functions/spatial_containers/point_object/PointObject<TObject>/UpdatePoint.md)
- [StlBvhTree](../../../../functions/spatial_containers/stl_bvh_tree/StlBvhTree.md)
- [EquationIdVector](../../../../functions/templates/classes/condition_template/EquationIdVector.md)
- [GetDofList](../../../../functions/templates/classes/condition_template/GetDofList.md)
- [Check](../../../../functions/templates/classes/condition_template/Check.md)
- [EquationIdVector](../../../../functions/templates/classes/element_template/EquationIdVector.md)
- [GetDofList](../../../../functions/templates/classes/element_template/GetDofList.md)
- [Check](../../../../functions/templates/classes/element_template/Check.md)
- [GetShape](../../../../functions/tensor_adaptors/connectivity_ids_tensor_adaptor/ConnectivityIdsTensorAdaptor/GetShape.md)
- [CheckContainer](../../../../functions/tensor_adaptors/connectivity_ids_tensor_adaptor/ConnectivityIdsTensorAdaptor/CheckContainer.md)
- [CollectIds](../../../../functions/tensor_adaptors/connectivity_ids_tensor_adaptor/ConnectivityIdsTensorAdaptor/CollectIds.md)
- [GetIntegrationMethod](../../../../functions/tensor_adaptors/geometries_tensor_adaptor/GeometriesTensorAdaptor/GetIntegrationMethod.md)
- [GetShape](../../../../functions/tensor_adaptors/geometries_tensor_adaptor/GeometriesTensorAdaptor/GetShape.md)
- [CollectIntegrationWeights](../../../../functions/tensor_adaptors/geometries_tensor_adaptor/GeometriesTensorAdaptor/CollectIntegrationWeights.md)
- [CollectShapeFunctions](../../../../functions/tensor_adaptors/geometries_tensor_adaptor/GeometriesTensorAdaptor/CollectShapeFunctions.md)
- [CollectShapeFunctionsDerivatives](../../../../functions/tensor_adaptors/geometries_tensor_adaptor/GeometriesTensorAdaptor/CollectShapeFunctionsDerivatives.md)
- [CollectNormals](../../../../functions/tensor_adaptors/geometries_tensor_adaptor/GeometriesTensorAdaptor/CollectNormals.md)
- [CollectJacobians](../../../../functions/tensor_adaptors/geometries_tensor_adaptor/GeometriesTensorAdaptor/CollectJacobians.md)
- [FillDomainSize](../../../../functions/tensor_adaptors/geometry_metrics_tensor_adaptor/FillDomainSize.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/containers/test_model_part_geometry_container/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/elements/test_edge_based_gradient_recovery_element/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/elements/test_levelset_convection_element_simplex/KRATOS_TEST_CASE_IN_SUITE.md)
- [validatewithtolerance](../../../../functions/tests/cpp_tests/external_libraries/test_delaunator/validatewithtolerance.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/geometries/test_coupling_geometry/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/geometries/test_nurbs_curve/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/geometries/test_nurbs_surface/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/geometries/test_nurbs_volume/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/modeler/test_voxel_mesh_modeler/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/processes/test_calculate_discontinuous_distance_to_skin_process/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/processes/test_coarsening_process/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/processes/test_find_nodal_neighbours_process/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/processes/test_mortar_mapper_process/KRATOS_TEST_CASE_IN_SUITE.md)
- [SetUpDistanceField](../../../../functions/tests/cpp_tests/processes/test_parallel_distance_calculation_process/SetUpDistanceField.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/processes/test_replace_elements_and_condition_process/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/processes/test_set_initial_state_process/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/processes/test_structured_mesh_generator_process/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/sources/test_condition/KRATOS_TEST_CASE_IN_SUITE.md)
- [EquationIdVector](../../../../functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/PrimalElement/EquationIdVector.md)
- [GetDofList](../../../../functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/PrimalElement/GetDofList.md)
- [GetValuesVector](../../../../functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/PrimalElement/GetValuesVector.md)
- [GetFirstDerivativesVector](../../../../functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/PrimalElement/GetFirstDerivativesVector.md)
- [GetSecondDerivativesVector](../../../../functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/PrimalElement/GetSecondDerivativesVector.md)
- [CalculateLeftHandSide](../../../../functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/PrimalElement/CalculateLeftHandSide.md)
- [CalculateRightHandSide](../../../../functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/PrimalElement/CalculateRightHandSide.md)
- [EquationIdVector](../../../../functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/PrimalCondition/EquationIdVector.md)
- [GetDofList](../../../../functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/PrimalCondition/GetDofList.md)
- [GetValuesVector](../../../../functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/PrimalCondition/GetValuesVector.md)
- [GetFirstDerivativesVector](../../../../functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/PrimalCondition/GetFirstDerivativesVector.md)
- [GetSecondDerivativesVector](../../../../functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/PrimalCondition/GetSecondDerivativesVector.md)
- [CalculateLeftHandSide](../../../../functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/PrimalCondition/CalculateLeftHandSide.md)
- [CalculateRightHandSide](../../../../functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/PrimalCondition/CalculateRightHandSide.md)
- [GetFirstDerivativesVector](../../../../functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/ThisExtensions/GetFirstDerivativesVector.md)
- [GetSecondDerivativesVector](../../../../functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/ThisExtensions/GetSecondDerivativesVector.md)
- [GetAuxiliaryVector](../../../../functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/ThisExtensions/GetAuxiliaryVector.md)
- [EquationIdVector](../../../../functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/AdjointElement/EquationIdVector.md)
- [GetDofList](../../../../functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/AdjointElement/GetDofList.md)
- [GetValuesVector](../../../../functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/AdjointElement/GetValuesVector.md)
- [GetFirstDerivativesVector](../../../../functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/AdjointElement/GetFirstDerivativesVector.md)
- [GetSecondDerivativesVector](../../../../functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/AdjointElement/GetSecondDerivativesVector.md)
- [CalculateSensitivityMatrix](../../../../functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/AdjointElement/CalculateSensitivityMatrix.md)
- [EquationIdVector](../../../../functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/AdjointCondition/EquationIdVector.md)
- [GetDofList](../../../../functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/AdjointCondition/GetDofList.md)
- [GetValuesVector](../../../../functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/AdjointCondition/GetValuesVector.md)
- [GetFirstDerivativesVector](../../../../functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/AdjointCondition/GetFirstDerivativesVector.md)
- [GetSecondDerivativesVector](../../../../functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/AdjointCondition/GetSecondDerivativesVector.md)
- [CalculateGradient](../../../../functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/ResponseFunction/CalculateGradient.md)
- [CalculateFirstDerivativesGradient](../../../../functions/tests/cpp_tests/strategies/schemes/test_residual_based_adjoint_bossak_scheme/ResponseFunction/CalculateFirstDerivativesGradient.md)
- [GetDofList](../../../../functions/tests/cpp_tests/strategies/strategies/explicit_strategy_testing_utils/AuxiliaryExplicitStrategiesTestElement/GetDofList.md)
- [AddExplicitContribution](../../../../functions/tests/cpp_tests/strategies/strategies/explicit_strategy_testing_utils/AuxiliaryExplicitStrategiesTestElement/AddExplicitContribution.md)
- [EquationIdVector](../../../../functions/tests/cpp_tests/strategies/strategies/explicit_strategy_testing_utils/AuxiliaryExplicitStrategiesTestElement/EquationIdVector.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/utilities/test_brute_force_point_locator/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/utilities/test_delaunator_utilities/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/utilities/test_discont_utils/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/utilities/test_divide_tetrahedra_3d_4/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/utilities/test_divide_triangle_2d_3/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/utilities/test_divide_triangle_3d_3/KRATOS_TEST_CASE_IN_SUITE.md)
- [MLSShapeFunctionsUtilityConsistencyCheck](../../../../functions/tests/cpp_tests/utilities/test_mls_shape_functions_utility/MLSShapeFunctionsUtilityConsistencyCheck.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/utilities/test_model_part_utils/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/utilities/test_sensitivity_builder/KRATOS_TEST_CASE_IN_SUITE.md)
- [test_cad_json_input_read](../../../../functions/tests/test_cad_json_input/TestCadJsonInput/test_cad_json_input_read.md)
- [SetupModelPart](../../../../functions/tests/test_coordinate_transformation_utils/SetupModelPart.md)
- [test_line_exact_integration_1](../../../../functions/tests/test_exact_integration/TestExactIntegration/test_line_exact_integration_1.md)
- [test_line_exact_integration_2](../../../../functions/tests/test_exact_integration/TestExactIntegration/test_line_exact_integration_2.md)
- [test_line_exact_integration_3](../../../../functions/tests/test_exact_integration/TestExactIntegration/test_line_exact_integration_3.md)
- [_test_triangle_exact_integration_1](../../../../functions/tests/test_exact_integration/TestExactIntegration/_test_triangle_exact_integration_1.md)
- [_test_triangle_exact_integration_2](../../../../functions/tests/test_exact_integration/TestExactIntegration/_test_triangle_exact_integration_2.md)
- [_test_triangle_exact_integration_3](../../../../functions/tests/test_exact_integration/TestExactIntegration/_test_triangle_exact_integration_3.md)
- [_test_quadrilateral_exact_integration_1](../../../../functions/tests/test_exact_integration/TestExactIntegration/_test_quadrilateral_exact_integration_1.md)
- [_test_quadrilateral_exact_integration_2](../../../../functions/tests/test_exact_integration/TestExactIntegration/_test_quadrilateral_exact_integration_2.md)
- [test_geometry_data](../../../../functions/tests/test_geometry_data/TestGeometryData/test_geometry_data.md)
- [test_geometry_container](../../../../functions/tests/test_model_part/TestModelPart/test_geometry_container.md)
- [test_Substract](../../../../functions/tests/test_model_part_operation_utilities/TestModelPartOperationUtilities/test_Substract.md)
- [test_ModelPartUtils_FromConnectivityGenerateElements](../../../../functions/tests/test_model_part_utils_connectivity_generations/TestModelPartUtilsConnectivityGenerations/test_ModelPartUtils_FromConnectivityGenerateElements.md)
- [test_ModelPartUtils_FromConnectivityGenerateConditions](../../../../functions/tests/test_model_part_utils_connectivity_generations/TestModelPartUtilsConnectivityGenerations/test_ModelPartUtils_FromConnectivityGenerateConditions.md)
- [__generate_discontinous_case](../../../../functions/tests/test_mortar_mapper/TestMortarMapperCore/__generate_discontinous_case.md)
- [FiniteDifferenceNormalShapeSensitivityTest](../../../../functions/tests/test_normal_utils/FiniteDifferenceNormalShapeSensitivityTest.md)
- [debug_vtk](../../../../functions/tests/test_obj_io/debug_vtk.md)
- [test_WriteObjIO](../../../../functions/tests/test_obj_io/TestObjIO/test_WriteObjIO.md)
- [test_assign_vector_variable_to_conditions_analytic_direction](../../../../functions/tests/test_processes/TestProcesses/test_assign_vector_variable_to_conditions_analytic_direction.md)
- [test_split_internal_interfaces_process_3D](../../../../functions/tests/test_processes/TestProcesses/test_split_internal_interfaces_process_3D.md)
- [test_parallel_redistance_maintain_plane_3d](../../../../functions/tests/test_redistance/TestRedistance/test_parallel_redistance_maintain_plane_3d.md)
- [CalculateEntityNeighbourMaps](../../../../functions/tests/test_sensitivity_utilities/CalculateEntityNeighbourMaps.md)
- [__InitializeAssignEntityDerivativesToNodesValues](../../../../functions/tests/test_sensitivity_utilities/TestSensitivityUtilitiesTwoDimSymmetricalSquare/__InitializeAssignEntityDerivativesToNodesValues.md)
- [test_WriteStlIO](../../../../functions/tests/test_stl_io/TestStlIO/test_WriteStlIO.md)
- [testGeometryMetricsTensorAdaptorDomainSizeCondition](../../../../functions/tests/test_tensor_adaptors/TestTensorAdaptors/testGeometryMetricsTensorAdaptorDomainSizeCondition.md)
- [testGeometryMetricsTensorAdaptorDomainSizeElement](../../../../functions/tests/test_tensor_adaptors/TestTensorAdaptors/testGeometryMetricsTensorAdaptorDomainSizeElement.md)
- [Create](../../../../functions/tests/test_utilities/test_bar_element/TestBarElement/Create.md)
- [Clone](../../../../functions/tests/test_utilities/test_bar_element/TestBarElement/Clone.md)
- [GetDofList](../../../../functions/tests/test_utilities/test_bar_element/TestBarElement/GetDofList.md)
- [EquationIdVector](../../../../functions/tests/test_utilities/test_bar_element/TestBarElement/EquationIdVector.md)
- [GetValuesVector](../../../../functions/tests/test_utilities/test_bar_element/TestBarElement/GetValuesVector.md)
- [GetFirstDerivativesVector](../../../../functions/tests/test_utilities/test_bar_element/TestBarElement/GetFirstDerivativesVector.md)
- [GetSecondDerivativesVector](../../../../functions/tests/test_utilities/test_bar_element/TestBarElement/GetSecondDerivativesVector.md)
- [CalculateLeftHandSide](../../../../functions/tests/test_utilities/test_bar_element/TestBarElement/CalculateLeftHandSide.md)
- [CalculateMassMatrix](../../../../functions/tests/test_utilities/test_bar_element/TestBarElement/CalculateMassMatrix.md)
- [CalculateLumpedMassVector](../../../../functions/tests/test_utilities/test_bar_element/TestBarElement/CalculateLumpedMassVector.md)
- [Check](../../../../functions/tests/test_utilities/test_bar_element/TestBarElement/Check.md)
- [Create](../../../../functions/tests/test_utilities/test_element/TestElement/Create.md)
- [Clone](../../../../functions/tests/test_utilities/test_element/TestElement/Clone.md)
- [GetDofList](../../../../functions/tests/test_utilities/test_element/TestElement/GetDofList.md)
- [EquationIdVector](../../../../functions/tests/test_utilities/test_element/TestElement/EquationIdVector.md)
- [GetValuesVector](../../../../functions/tests/test_utilities/test_element/TestElement/GetValuesVector.md)
- [GetFirstDerivativesVector](../../../../functions/tests/test_utilities/test_element/TestElement/GetFirstDerivativesVector.md)
- [GetSecondDerivativesVector](../../../../functions/tests/test_utilities/test_element/TestElement/GetSecondDerivativesVector.md)
- [CalculateRightHandSide](../../../../functions/tests/test_utilities/test_element/TestElement/CalculateRightHandSide.md)
- [CalculateLeftHandSide](../../../../functions/tests/test_utilities/test_element/TestElement/CalculateLeftHandSide.md)
- [CalculateMassMatrix](../../../../functions/tests/test_utilities/test_element/TestElement/CalculateMassMatrix.md)
- [CalculateDampingMatrix](../../../../functions/tests/test_utilities/test_element/TestElement/CalculateDampingMatrix.md)
- [Check](../../../../functions/tests/test_utilities/test_element/TestElement/Check.md)
- [Initialize](../../../../functions/tests/test_utilities/test_element/TestElement/Initialize.md)
- [Create](../../../../functions/tests/test_utilities/test_laplacian_element/TestLaplacianElement/Create.md)
- [Clone](../../../../functions/tests/test_utilities/test_laplacian_element/TestLaplacianElement/Clone.md)
- [GetDofList](../../../../functions/tests/test_utilities/test_laplacian_element/TestLaplacianElement/GetDofList.md)
- [EquationIdVector](../../../../functions/tests/test_utilities/test_laplacian_element/TestLaplacianElement/EquationIdVector.md)
- [GetValuesVector](../../../../functions/tests/test_utilities/test_laplacian_element/TestLaplacianElement/GetValuesVector.md)
- [CalculateLocalSystem](../../../../functions/tests/test_utilities/test_laplacian_element/TestLaplacianElement/CalculateLocalSystem.md)
- [CalculateMassMatrix](../../../../functions/tests/test_utilities/test_laplacian_element/TestLaplacianElement/CalculateMassMatrix.md)
- [CalculateDampingMatrix](../../../../functions/tests/test_utilities/test_laplacian_element/TestLaplacianElement/CalculateDampingMatrix.md)
- [Check](../../../../functions/tests/test_utilities/test_laplacian_element/TestLaplacianElement/Check.md)
- [test_nodal_area](../../../../functions/tests/test_variable_redistribution/TestVariableRedistributionUtility/test_nodal_area.md)
- [test_TestVoxelMeshModelerNoMdpa](../../../../functions/tests/test_voxel_modeler/TestVoxelMeshModeler/test_TestVoxelMeshModelerNoMdpa.md)
- [test_TestVoxelMeshModelerElementInput](../../../../functions/tests/test_voxel_modeler/TestVoxelMeshModeler/test_TestVoxelMeshModelerElementInput.md)
- [test_TestVoxelMeshModelerQuadraticElementInput](../../../../functions/tests/test_voxel_modeler/TestVoxelMeshModeler/test_TestVoxelMeshModelerQuadraticElementInput.md)
- [test_TestVoxelMeshModelerElementInputAutoSkin](../../../../functions/tests/test_voxel_modeler/TestVoxelMeshModeler/test_TestVoxelMeshModelerElementInputAutoSkin.md)
- [test_TestVoxelMeshModelerTetrahedralElementInput](../../../../functions/tests/test_voxel_modeler/TestVoxelMeshModeler/test_TestVoxelMeshModelerTetrahedralElementInput.md)
- [test_TestVoxelMeshModelerInRadius](../../../../functions/tests/test_voxel_modeler/TestVoxelMeshModeler/test_TestVoxelMeshModelerInRadius.md)
- [test_TestWallConditionVoxelizer](../../../../functions/tests/test_voxel_modeler/TestBoundaryConditionVoxelizer/test_TestWallConditionVoxelizer.md)
- [test_TestFacesBetweenColors](../../../../functions/tests/test_voxel_modeler/TestBoundaryConditionVoxelizer/test_TestFacesBetweenColors.md)
- [test_TestOuterWallConditionVoxelizer](../../../../functions/tests/test_voxel_modeler/TestBoundaryConditionVoxelizer/test_TestOuterWallConditionVoxelizer.md)
- [fill_nodes](../../../../functions/tests/test_vtu_output/fill_nodes.md)
- [GetNodes](../../../../functions/tests/test_vtu_output/TestVtuOutput/GetNodes.md)
- [GetNumberOfNodes](../../../../functions/tests/test_vtu_output/TestVtuOutput/GetNumberOfNodes.md)
- [CheckGaussVtuFile](../../../../functions/tests/test_vtu_output/TestVtuOutput/CheckGaussVtuFile.md)
- [ActivateElementsAndConditions](../../../../functions/utilities/activation_utilities/ActivationUtilities/ActivateElementsAndConditions.md)
- [AddElementWithNodes](../../../../functions/utilities/auxiliar_model_part_utilities/AuxiliarModelPartUtilities/AddElementWithNodes.md)
- [AddConditionWithNodes](../../../../functions/utilities/auxiliar_model_part_utilities/AuxiliarModelPartUtilities/AddConditionWithNodes.md)
- [RemoveElementAndBelongings](../../../../functions/utilities/auxiliar_model_part_utilities/AuxiliarModelPartUtilities/RemoveElementAndBelongings.md)
- [RemoveElementsAndBelongings](../../../../functions/utilities/auxiliar_model_part_utilities/AuxiliarModelPartUtilities/RemoveElementsAndBelongings.md)
- [RemoveConditionAndBelongings](../../../../functions/utilities/auxiliar_model_part_utilities/AuxiliarModelPartUtilities/RemoveConditionAndBelongings.md)
- [RemoveConditionsAndBelongings](../../../../functions/utilities/auxiliar_model_part_utilities/AuxiliarModelPartUtilities/RemoveConditionsAndBelongings.md)
- [RemoveOrphanNodesFromSubModelParts](../../../../functions/utilities/auxiliar_model_part_utilities/AuxiliarModelPartUtilities/RemoveOrphanNodesFromSubModelParts.md)
- [RetrieveElementsNeighbourElementsIds](../../../../functions/utilities/auxiliar_model_part_utilities/AuxiliarModelPartUtilities/RetrieveElementsNeighbourElementsIds.md)
- [RetrieveConditionsNeighbourConditionsIds](../../../../functions/utilities/auxiliar_model_part_utilities/AuxiliarModelPartUtilities/RetrieveConditionsNeighbourConditionsIds.md)
- [AuxiliaryAddEntitiesWithNodes](../../../../functions/utilities/auxiliar_model_part_utilities/AuxiliaryAddEntitiesWithNodes.md)
- [FindPointOnMesh](../../../../functions/utilities/binbased_fast_point_locator/FindPointOnMesh.md)
- [FindNodesInElement](../../../../functions/utilities/binbased_nodes_in_element_locator/BinBasedNodesInElementLocator/FindNodesInElement.md)
- [CalculateDistances](../../../../functions/utilities/body_distance_calculation_utils/CalculateDistances.md)
- [CalculateBodyNormals](../../../../functions/utilities/body_normal_calculation_utils/BodyNormalCalculationUtils/CalculateBodyNormals.md)
- [FindObject](../../../../functions/utilities/brute_force_point_locator/BruteForcePointLocator/FindObject.md)
- [MoveParticles_Substepping](../../../../functions/utilities/convect_particles_utilities/ParticleConvectUtily/MoveParticles_Substepping.md)
- [MoveParticles_RK4](../../../../functions/utilities/convect_particles_utilities/ParticleConvectUtily/MoveParticles_RK4.md)
- [EraseOuterElements](../../../../functions/utilities/convect_particles_utilities/ParticleConvectUtily/EraseOuterElements.md)
- [FindSmallestEdge](../../../../functions/utilities/cutting_utility/CuttingUtility/FindSmallestEdge.md)
- [AddSkinConditions](../../../../functions/utilities/cutting_utility/CuttingUtility/AddSkinConditions.md)
- [FirstLoop](../../../../functions/utilities/cutting_utility/CuttingUtility/FirstLoop.md)
- [GenerateElements](../../../../functions/utilities/cutting_utility/CuttingUtility/GenerateElements.md)
- [GenerateSkin](../../../../functions/utilities/embedded_skin_utility/EmbeddedSkinUtility<TDim>/GenerateSkin.md)
- [ComputeElementSkin](../../../../functions/utilities/embedded_skin_utility/EmbeddedSkinUtility<TDim>/ComputeElementSkin.md)
- [SetDistancesVector](../../../../functions/utilities/embedded_skin_utility/EmbeddedSkinUtility<TDim>/SetDistancesVector.md)
- [InterpolateMeshVariableToSkinSpecialization](../../../../functions/utilities/embedded_skin_utility/InterpolateMeshVariableToSkinSpecialization.md)
- [GenerateSingleType](../../../../functions/utilities/entities_utilities/EntitityIdentifier<TEntity>/GenerateSingleType.md)
- [GenerateMultipleTypes](../../../../functions/utilities/entities_utilities/EntitityIdentifier<TEntity>/GenerateMultipleTypes.md)
- [GenerateTemplatedTypes](../../../../functions/utilities/entities_utilities/EntitityIdentifier<TEntity>/GenerateTemplatedTypes.md)
- [TestGetExactIntegration](../../../../functions/utilities/exact_mortar_segmentation_utility/ExactMortarIntegrationUtility<TDim, TNumNodes, TBelong, TNumNodesMaster>/TestGetExactIntegration.md)
- [TestGetExactAreaIntegration](../../../../functions/utilities/exact_mortar_segmentation_utility/ExactMortarIntegrationUtility<TDim, TNumNodes, TBelong, TNumNodesMaster>/TestGetExactAreaIntegration.md)
- [TestIODebug](../../../../functions/utilities/exact_mortar_segmentation_utility/ExactMortarIntegrationUtility<TDim, TNumNodes, TBelong, TNumNodesMaster>/TestIODebug.md)
- [SkippableEntity](../../../../functions/utilities/input_output_utilities/InputOutputUtilities/SkippableEntity.md)
- [AddSkinConditions](../../../../functions/utilities/iso_printer/IsosurfacePrinterApplication/AddSkinConditions.md)
- [AddScalarVarIsosurface](../../../../functions/utilities/iso_printer/IsosurfacePrinterApplication/AddScalarVarIsosurface.md)
- [AddScalarVarIsosurfaceAndLower](../../../../functions/utilities/iso_printer/IsosurfacePrinterApplication/AddScalarVarIsosurfaceAndLower.md)
- [AddScalarVarIsosurfaceAndHigher](../../../../functions/utilities/iso_printer/IsosurfacePrinterApplication/AddScalarVarIsosurfaceAndHigher.md)
- [ComputeGraph](../../../../functions/utilities/model_part_graph_utilities/ModelPartGraphUtilities/ComputeGraph.md)
- [FillNodePointersForEntities](../../../../functions/utilities/model_part_operation_utilities/ModelPartOperationUtilities/FillNodePointersForEntities.md)
- [FillNodesFromEntities](../../../../functions/utilities/model_part_operation_utilities/ModelPartOperationUtilities/FillNodesFromEntities.md)
- [FindNeighbourEntities](../../../../functions/utilities/model_part_operation_utilities/ModelPartOperationUtilities/FindNeighbourEntities.md)
- [CheckEntities](../../../../functions/utilities/model_part_operation_utilities/ModelPartOperationUtilities/CheckEntities.md)
- [AddEntities](../../../../functions/utilities/model_part_operation_utilities/AddEntities.md)
- [GenerateEntitiesFromConnectivities](../../../../functions/utilities/model_part_utils/GenerateEntitiesFromConnectivities.md)
- [InvertNormalForFlag](../../../../functions/utilities/mortar_utilities/InvertNormalForFlag.md)
- [CalculateNormalsInContainer](../../../../functions/utilities/normal_calculation_utils/NormalCalculationUtils/CalculateNormalsInContainer.md)
- [InitializeNormals](../../../../functions/utilities/normal_calculation_utils/NormalCalculationUtils/InitializeNormals.md)
- [CalculateNormalShapeDerivativesOnSimplex](../../../../functions/utilities/normal_calculation_utils/NormalCalculationUtils/CalculateNormalShapeDerivativesOnSimplex.md)
- [for](../../../../functions/utilities/normal_calculation_utils/for.md)
- [ComputeUnitNormalsFromAreaNormals](../../../../functions/utilities/normal_calculation_utils/NormalCalculationUtils/ComputeUnitNormalsFromAreaNormals.md)
- [CalculateNormal2D](../../../../functions/utilities/normal_calculation_utils/NormalCalculationUtils/CalculateNormal2D.md)
- [CalculateNormal3D](../../../../functions/utilities/normal_calculation_utils/NormalCalculationUtils/CalculateNormal3D.md)
- [CalculateNormalShapeDerivative3D](../../../../functions/utilities/normal_calculation_utils/NormalCalculationUtils/CalculateNormalShapeDerivative3D.md)
- [CalculateNormalsUsingGenericAlgorithm](../../../../functions/utilities/normal_calculation_utils/NormalCalculationUtils/CalculateNormalsUsingGenericAlgorithm.md)
- [CheckUseSimplex](../../../../functions/utilities/normal_calculation_utils/NormalCalculationUtils/CheckUseSimplex.md)
- [AuxiliaryCalculateOnSimplex](../../../../functions/utilities/normal_calculation_utils/NormalCalculationUtils/AuxiliaryCalculateOnSimplex.md)
- [CalculateOnSimplex](../../../../functions/utilities/normal_calculation_utils/CalculateOnSimplex.md)
- [CalculateOnSimplexLowMemory](../../../../functions/utilities/normal_calculation_utils/CalculateOnSimplexLowMemory.md)
- [CountParticlesInNodes](../../../../functions/utilities/particles_utilities/CountParticlesInNodes.md)
- [InterpolateValuesAtCoordinates](../../../../functions/utilities/particles_utilities/InterpolateValuesAtCoordinates.md)
- [QuadraticErrorFunctionPoint](../../../../functions/utilities/qef_utility/QuadraticErrorFunction/QuadraticErrorFunctionPoint.md)
- [AuxiliaryAssignPeriodicity](../../../../functions/utilities/rve_periodicity_utility/RVEPeriodicityUtility/AuxiliaryAssignPeriodicity.md)
- [AssembleContainerContributions](../../../../functions/utilities/sensitivity_builder/AssembleContainerContributions.md)
- [AssignEntityDerivativesToNodes](../../../../functions/utilities/sensitivity_utilities/SensitivityUtilities/AssignEntityDerivativesToNodes.md)
- [CalculateGradientBasedConformingExtensionBasis](../../../../functions/utilities/shifted_boundary_meshless_interface_utility/ShiftedBoundaryMeshlessInterfaceUtility/CalculateGradientBasedConformingExtensionBasis.md)
- [CalculateMeshlessBasedConformingExtensionBasis](../../../../functions/utilities/shifted_boundary_meshless_interface_utility/ShiftedBoundaryMeshlessInterfaceUtility/CalculateMeshlessBasedConformingExtensionBasis.md)
- [CalculateMeshlessBasedNonConformingExtensionBasis](../../../../functions/utilities/shifted_boundary_meshless_interface_utility/ShiftedBoundaryMeshlessInterfaceUtility/CalculateMeshlessBasedNonConformingExtensionBasis.md)
- [SetInterfaceFlags](../../../../functions/utilities/shifted_boundary_meshless_interface_utility/ShiftedBoundaryMeshlessInterfaceUtility/SetInterfaceFlags.md)
- [SetSplitElementSupportCloud](../../../../functions/utilities/shifted_boundary_meshless_interface_utility/ShiftedBoundaryMeshlessInterfaceUtility/SetSplitElementSupportCloud.md)
- [SetSurrogateBoundaryNodalGradientWeights](../../../../functions/utilities/shifted_boundary_meshless_interface_utility/ShiftedBoundaryMeshlessInterfaceUtility/SetSurrogateBoundaryNodalGradientWeights.md)
- [CalculateDistances2D](../../../../functions/utilities/signed_distance_calculation_utils/SignedDistanceCalculationUtils/CalculateDistances2D.md)
- [CalculateDistances3D](../../../../functions/utilities/signed_distance_calculation_utils/SignedDistanceCalculationUtils/CalculateDistances3D.md)
- [CalculateBoundingBox](../../../../functions/utilities/spatial_containers_configure/SpatialContainersConfigure/CalculateBoundingBox.md)
- [Intersection](../../../../functions/utilities/spatial_containers_configure/SpatialContainersConfigure/Intersection.md)
- [IntersectionBox](../../../../functions/utilities/spatial_containers_configure/SpatialContainersConfigure/IntersectionBox.md)
- [DetermineIfCompatibleGeometries](../../../../functions/utilities/specifications_utilities/SpecificationsUtilities/DetermineIfCompatibleGeometries.md)
- [CreateTriangleMeshFromNodes](../../../../functions/utilities/tessellation_utilities/delaunator_utilities/CreateTriangleMeshFromNodes.md)
- [CallSpecializedConvertDistributedValuesToPoint](../../../../functions/utilities/variable_redistribution_utility/VariableRedistributionUtility/CallSpecializedConvertDistributedValuesToPoint.md)
- [CallSpecializedDistributePointValues](../../../../functions/utilities/variable_redistribution_utility/VariableRedistributionUtility/CallSpecializedDistributePointValues.md)
- [SpecializedConvertDistributedValuesToPoint](../../../../functions/utilities/variable_redistribution_utility/VariableRedistributionUtility/SpecializedConvertDistributedValuesToPoint.md)
- [ComputeNodalSizes](../../../../functions/utilities/variable_redistribution_utility/VariableRedistributionUtility/ComputeNodalSizes.md)
- [UpdateDistributionRHS](../../../../functions/utilities/variable_redistribution_utility/VariableRedistributionUtility/UpdateDistributionRHS.md)
- [WeightedAccumulateVariableOnNodes](../../../../functions/utilities/variable_utils/VariableUtils/WeightedAccumulateVariableOnNodes.md)