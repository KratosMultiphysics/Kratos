---
type: C++ Method
title: CreateNewProperties
resource: sources/model_part.cpp#L670-L689
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/model_part/IsSubModelPart
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/input_output/obj_io/ObjIO/ReadModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/stl_io/StlIO/ReadModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/entity_generation/voxel_mesher_entity_generation/VoxelMesherEntityGeneration/CreateAndGetProperty
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/tests/cpp_tests/sources/test_mpi_communicator/ModelPartForMPICommunicatorTests
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/tests/test_gather_modelpart_utility/CreateLocalModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/tests/test_utilities/mpi_cpp_test_utilities/MPICppTestUtilities/GenerateDistributedBarStructure
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/tests/test_utilities/mpi_cpp_test_utilities/MPICppTestUtilities/GenerateDistributedTriangleMesh
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/tests/test_utilities/mpi_cpp_test_utilities/MPICppTestUtilities/GenerateDistributedTriforceMesh
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_embedded_nodal_variable_from_skin_process/CalculateEmbeddedNodalVariableFromSkinProcess/GenerateIntersectedEdgesElementsModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/edge_based_gradient_recovery_process/EdgeBasedGradientRecoveryProcess/InitializeGradientRecoveryModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/find_intersected_geometrical_objects_with_obb_process/FindIntersectedGeometricalObjectsWithOBBProcess/FindIntersectedGeometricalObjectsWithOBBProcess
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/skin_detection_process/SkinDetectionProcess<TDim>/CreateConditions
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/split_internal_interfaces_process/SplitInternalInterfacesProcess/SplitBoundary
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/structured_mesh_generator_process/StructuredMeshGeneratorProcess/GenerateTriangularElements
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/structured_mesh_generator_process/StructuredMeshGeneratorProcess/GenerateTetrahedraElements
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/sub_model_part_skin_detection_process/SubModelPartSkinDetectionProcess<TDim>/CreateConditions
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python/add_model_part_to_python/AddModelPartToPython
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/containers/test_pointer_vector_set/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/containers/test_properties/KRATOS_TEST_CASE_IN_SUITE
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
  - target: functions/tests/cpp_tests/geometries/test_brep_surface/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/geometries/test_brep_volume/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/geometries/test_coupling_geometry/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/includes/test_property_accessor/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/modeler/test_copy_properties_modeler/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/processes/test_apply_ray_casting_interface_recognition_process/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/processes/test_calculate_discontinuous_distance_to_skin_process/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/processes/test_calculate_distance_to_skin_process/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/processes/test_fast_transfer_between_model_parts_process/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/processes/test_mortar_mapper_process/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/processes/test_replace_elements_and_condition_process/GenerateMixedModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/processes/test_replace_elements_and_condition_process/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/processes/test_set_initial_state_process/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/processes/test_skin_detection_process/CreateSimpleGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/sources/test_condition/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/sources/test_element/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/sources/test_model_part/GenerateGenericModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/sources/test_model_part/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/strategies/builder_and_solvers/test_builder_and_solver/BasicTestBuilderAndSolverDisplacement
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/strategies/builder_and_solvers/test_builder_and_solver/BasicTestBuilderAndSolverDisplacementAllDoFsMaster
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/strategies/builder_and_solvers/test_builder_and_solver/BasicTestBuilderAndSolverDisplacementWithZeroContribution
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/strategies/builder_and_solvers/test_builder_and_solver/ExtendedTestBuilderAndSolverDisplacement
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/strategies/builder_and_solvers/test_explicit_builder_and_solver/GenerateTestExplicitBuilderModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/strategies/strategies/test_strategies/BasicTestStrategyDisplacement
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/utilities/test_auxiliar_model_part_utilities/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/utilities/test_model_part_utils/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_connectivity_ids_tensor_adaptor/TestConnectivityIdsTensorAdaptor/setUpClass
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_connectivity_preserve_modeler/TestConnectivityPreserveModeler/test_properties_shared
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_dofs/TestDofs/testBlockBuildDofSetUtility
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_duplicate_mesh_modeler/TestDuplicateMeshModeler/test_DuplicateMeshModeler
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_geometries_tensor_adaptor/TestGeometriesTensorAdaptor/setUp
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_model_part/TestModelPart/test_model_part_sub_properties
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_model_part/TestModelPart/test_model_part_properties_container
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_model_part/TestModelPart/test_create_non_existing_element
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_model_part/TestModelPart/test_create_non_existing_condition
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_properties/TestProperties/test_copy_properties
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_properties/TestProperties/test_erase_properties
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_properties/TestProperties/test_properties_utilities
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_properties/TestProperties/test_has_table
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_properties/TestProperties/test_accessor
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_properties/TestProperties/test_accessor_invalid
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_properties/TestProperties/test_accessor_get_consumed
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_properties/TestProperties/test_accessor_assign_twice
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_properties/TestProperties/test_print
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_sparse_matrices/TestSparseMatrixInterface/test_get_equation_id_csr_indices_1
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_sparse_matrices/TestSparseMatrixInterface/test_assemble_with_csr_indices
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_sparse_matrices/TestSparseMatrixInterface/test_assemble_with_equation_ids
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_sparse_matrices/TestSparseMatrixInterface/test_assemble_with_csr_indices_int64
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_sparse_matrices/TestSparseMatrixInterface/test_assemble_with_equation_ids_int64
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_specifications_utilities/TestSpecificationsUtilities/test_specifications_utilities_GetDofsListFromSpecifications
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_tensor_adaptors/TestTensorAdaptors/setUpClass
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_utilities/cpp_tests_utilities/Create2DGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_utilities/cpp_tests_utilities/CreateTestModelPartTriangle2D3N
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_utilities/cpp_tests_utilities/Create2DQuadrilateralsGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_utilities/cpp_tests_utilities/Create3DGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_utilities/cpp_tests_utilities/CreateTestModelPartTetrahedra3D4N
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_utilities/cpp_tests_utilities/Create3DHexahedraGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_utilities/cpp_tests_utilities/Create3DQuadraticGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_vtu_output/TestVtuOutput/setUpClass
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/exact_mortar_segmentation_utility/ExactMortarIntegrationUtility<TDim, TNumNodes, TBelong, TNumNodesMaster>/TestIODebug
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/read_materials_utility/ReadMaterialsUtility/GetPropertyBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/read_materials_utility/ReadMaterialsUtility/CreateSubProperties
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/read_materials_utility/ReadMaterialsUtility/AssignPropertyBlock
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/tessellation_utilities/delaunator_utilities/CreateTriangleMeshFromNodes
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`ModelPart::PropertiesType::Pointer ModelPart::CreateNewProperties( IndexType PropertiesId, IndexType MeshIndex )`

# Calls

- [IsSubModelPart](../../../../functions/includes/model_part/IsSubModelPart.md)

# Called by

- [ReadModelPart](../../../../functions/input_output/obj_io/ObjIO/ReadModelPart.md)
- [ReadModelPart](../../../../functions/input_output/stl_io/StlIO/ReadModelPart.md)
- [CreateAndGetProperty](../../../../functions/modeler/entity_generation/voxel_mesher_entity_generation/VoxelMesherEntityGeneration/CreateAndGetProperty.md)
- [ModelPartForMPICommunicatorTests](../../../../functions/mpi/tests/cpp_tests/sources/test_mpi_communicator/ModelPartForMPICommunicatorTests.md)
- [CreateLocalModelPart](../../../../functions/mpi/tests/test_gather_modelpart_utility/CreateLocalModelPart.md)
- [GenerateDistributedBarStructure](../../../../functions/mpi/tests/test_utilities/mpi_cpp_test_utilities/MPICppTestUtilities/GenerateDistributedBarStructure.md)
- [GenerateDistributedTriangleMesh](../../../../functions/mpi/tests/test_utilities/mpi_cpp_test_utilities/MPICppTestUtilities/GenerateDistributedTriangleMesh.md)
- [GenerateDistributedTriforceMesh](../../../../functions/mpi/tests/test_utilities/mpi_cpp_test_utilities/MPICppTestUtilities/GenerateDistributedTriforceMesh.md)
- [GenerateIntersectedEdgesElementsModelPart](../../../../functions/processes/calculate_embedded_nodal_variable_from_skin_process/CalculateEmbeddedNodalVariableFromSkinProcess/GenerateIntersectedEdgesElementsModelPart.md)
- [InitializeGradientRecoveryModelPart](../../../../functions/processes/edge_based_gradient_recovery_process/EdgeBasedGradientRecoveryProcess/InitializeGradientRecoveryModelPart.md)
- [FindIntersectedGeometricalObjectsWithOBBProcess](../../../../functions/processes/find_intersected_geometrical_objects_with_obb_process/FindIntersectedGeometricalObjectsWithOBBProcess/FindIntersectedGeometricalObjectsWithOBBProcess.md)
- [CreateConditions](../../../../functions/processes/skin_detection_process/SkinDetectionProcess<TDim>/CreateConditions.md)
- [SplitBoundary](../../../../functions/processes/split_internal_interfaces_process/SplitInternalInterfacesProcess/SplitBoundary.md)
- [GenerateTriangularElements](../../../../functions/processes/structured_mesh_generator_process/StructuredMeshGeneratorProcess/GenerateTriangularElements.md)
- [GenerateTetrahedraElements](../../../../functions/processes/structured_mesh_generator_process/StructuredMeshGeneratorProcess/GenerateTetrahedraElements.md)
- [CreateConditions](../../../../functions/processes/sub_model_part_skin_detection_process/SubModelPartSkinDetectionProcess<TDim>/CreateConditions.md)
- [AddModelPartToPython](../../../../functions/python/add_model_part_to_python/AddModelPartToPython.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/containers/test_pointer_vector_set/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/containers/test_properties/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/elements/test_edge_based_gradient_recovery_element/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/elements/test_levelset_convection_element_simplex/KRATOS_TEST_CASE_IN_SUITE.md)
- [validatewithtolerance](../../../../functions/tests/cpp_tests/external_libraries/test_delaunator/validatewithtolerance.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/geometries/test_brep_surface/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/geometries/test_brep_volume/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/geometries/test_coupling_geometry/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/includes/test_property_accessor/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/modeler/test_copy_properties_modeler/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/processes/test_apply_ray_casting_interface_recognition_process/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/processes/test_calculate_discontinuous_distance_to_skin_process/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/processes/test_calculate_distance_to_skin_process/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/processes/test_fast_transfer_between_model_parts_process/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/processes/test_mortar_mapper_process/KRATOS_TEST_CASE_IN_SUITE.md)
- [GenerateMixedModelPart](../../../../functions/tests/cpp_tests/processes/test_replace_elements_and_condition_process/GenerateMixedModelPart.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/processes/test_replace_elements_and_condition_process/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/processes/test_set_initial_state_process/KRATOS_TEST_CASE_IN_SUITE.md)
- [CreateSimpleGeometry](../../../../functions/tests/cpp_tests/processes/test_skin_detection_process/CreateSimpleGeometry.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/sources/test_condition/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/sources/test_element/KRATOS_TEST_CASE_IN_SUITE.md)
- [GenerateGenericModelPart](../../../../functions/tests/cpp_tests/sources/test_model_part/GenerateGenericModelPart.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/sources/test_model_part/KRATOS_TEST_CASE_IN_SUITE.md)
- [BasicTestBuilderAndSolverDisplacement](../../../../functions/tests/cpp_tests/strategies/builder_and_solvers/test_builder_and_solver/BasicTestBuilderAndSolverDisplacement.md)
- [BasicTestBuilderAndSolverDisplacementAllDoFsMaster](../../../../functions/tests/cpp_tests/strategies/builder_and_solvers/test_builder_and_solver/BasicTestBuilderAndSolverDisplacementAllDoFsMaster.md)
- [BasicTestBuilderAndSolverDisplacementWithZeroContribution](../../../../functions/tests/cpp_tests/strategies/builder_and_solvers/test_builder_and_solver/BasicTestBuilderAndSolverDisplacementWithZeroContribution.md)
- [ExtendedTestBuilderAndSolverDisplacement](../../../../functions/tests/cpp_tests/strategies/builder_and_solvers/test_builder_and_solver/ExtendedTestBuilderAndSolverDisplacement.md)
- [GenerateTestExplicitBuilderModelPart](../../../../functions/tests/cpp_tests/strategies/builder_and_solvers/test_explicit_builder_and_solver/GenerateTestExplicitBuilderModelPart.md)
- [BasicTestStrategyDisplacement](../../../../functions/tests/cpp_tests/strategies/strategies/test_strategies/BasicTestStrategyDisplacement.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/utilities/test_auxiliar_model_part_utilities/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/utilities/test_model_part_utils/KRATOS_TEST_CASE_IN_SUITE.md)
- [setUpClass](../../../../functions/tests/test_connectivity_ids_tensor_adaptor/TestConnectivityIdsTensorAdaptor/setUpClass.md)
- [test_properties_shared](../../../../functions/tests/test_connectivity_preserve_modeler/TestConnectivityPreserveModeler/test_properties_shared.md)
- [testBlockBuildDofSetUtility](../../../../functions/tests/test_dofs/TestDofs/testBlockBuildDofSetUtility.md)
- [test_DuplicateMeshModeler](../../../../functions/tests/test_duplicate_mesh_modeler/TestDuplicateMeshModeler/test_DuplicateMeshModeler.md)
- [setUp](../../../../functions/tests/test_geometries_tensor_adaptor/TestGeometriesTensorAdaptor/setUp.md)
- [test_model_part_sub_properties](../../../../functions/tests/test_model_part/TestModelPart/test_model_part_sub_properties.md)
- [test_model_part_properties_container](../../../../functions/tests/test_model_part/TestModelPart/test_model_part_properties_container.md)
- [test_create_non_existing_element](../../../../functions/tests/test_model_part/TestModelPart/test_create_non_existing_element.md)
- [test_create_non_existing_condition](../../../../functions/tests/test_model_part/TestModelPart/test_create_non_existing_condition.md)
- [test_copy_properties](../../../../functions/tests/test_properties/TestProperties/test_copy_properties.md)
- [test_erase_properties](../../../../functions/tests/test_properties/TestProperties/test_erase_properties.md)
- [test_properties_utilities](../../../../functions/tests/test_properties/TestProperties/test_properties_utilities.md)
- [test_has_table](../../../../functions/tests/test_properties/TestProperties/test_has_table.md)
- [test_accessor](../../../../functions/tests/test_properties/TestProperties/test_accessor.md)
- [test_accessor_invalid](../../../../functions/tests/test_properties/TestProperties/test_accessor_invalid.md)
- [test_accessor_get_consumed](../../../../functions/tests/test_properties/TestProperties/test_accessor_get_consumed.md)
- [test_accessor_assign_twice](../../../../functions/tests/test_properties/TestProperties/test_accessor_assign_twice.md)
- [test_print](../../../../functions/tests/test_properties/TestProperties/test_print.md)
- [test_get_equation_id_csr_indices_1](../../../../functions/tests/test_sparse_matrices/TestSparseMatrixInterface/test_get_equation_id_csr_indices_1.md)
- [test_assemble_with_csr_indices](../../../../functions/tests/test_sparse_matrices/TestSparseMatrixInterface/test_assemble_with_csr_indices.md)
- [test_assemble_with_equation_ids](../../../../functions/tests/test_sparse_matrices/TestSparseMatrixInterface/test_assemble_with_equation_ids.md)
- [test_assemble_with_csr_indices_int64](../../../../functions/tests/test_sparse_matrices/TestSparseMatrixInterface/test_assemble_with_csr_indices_int64.md)
- [test_assemble_with_equation_ids_int64](../../../../functions/tests/test_sparse_matrices/TestSparseMatrixInterface/test_assemble_with_equation_ids_int64.md)
- [test_specifications_utilities_GetDofsListFromSpecifications](../../../../functions/tests/test_specifications_utilities/TestSpecificationsUtilities/test_specifications_utilities_GetDofsListFromSpecifications.md)
- [setUpClass](../../../../functions/tests/test_tensor_adaptors/TestTensorAdaptors/setUpClass.md)
- [Create2DGeometry](../../../../functions/tests/test_utilities/cpp_tests_utilities/Create2DGeometry.md)
- [CreateTestModelPartTriangle2D3N](../../../../functions/tests/test_utilities/cpp_tests_utilities/CreateTestModelPartTriangle2D3N.md)
- [Create2DQuadrilateralsGeometry](../../../../functions/tests/test_utilities/cpp_tests_utilities/Create2DQuadrilateralsGeometry.md)
- [Create3DGeometry](../../../../functions/tests/test_utilities/cpp_tests_utilities/Create3DGeometry.md)
- [CreateTestModelPartTetrahedra3D4N](../../../../functions/tests/test_utilities/cpp_tests_utilities/CreateTestModelPartTetrahedra3D4N.md)
- [Create3DHexahedraGeometry](../../../../functions/tests/test_utilities/cpp_tests_utilities/Create3DHexahedraGeometry.md)
- [Create3DQuadraticGeometry](../../../../functions/tests/test_utilities/cpp_tests_utilities/Create3DQuadraticGeometry.md)
- [setUpClass](../../../../functions/tests/test_vtu_output/TestVtuOutput/setUpClass.md)
- [TestIODebug](../../../../functions/utilities/exact_mortar_segmentation_utility/ExactMortarIntegrationUtility<TDim, TNumNodes, TBelong, TNumNodesMaster>/TestIODebug.md)
- [GetPropertyBlock](../../../../functions/utilities/read_materials_utility/ReadMaterialsUtility/GetPropertyBlock.md)
- [CreateSubProperties](../../../../functions/utilities/read_materials_utility/ReadMaterialsUtility/CreateSubProperties.md)
- [AssignPropertyBlock](../../../../functions/utilities/read_materials_utility/ReadMaterialsUtility/AssignPropertyBlock.md)
- [CreateTriangleMeshFromNodes](../../../../functions/utilities/tessellation_utilities/delaunator_utilities/CreateTriangleMeshFromNodes.md)