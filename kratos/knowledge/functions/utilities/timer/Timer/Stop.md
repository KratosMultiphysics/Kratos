---
type: C++ Method
title: Stop
resource: utilities/timer.cpp#L125-L141
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/timer/CreateFullLabel
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/pointer_vector_set/PointerVectorSet/pop_back
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/timer/TimerData/GetStartTime
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/future/solving_strategies/schemes/implicit_scheme/ImplicitScheme/Build
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/future/solving_strategies/schemes/implicit_scheme/ImplicitScheme/BuildWithSafeAssemble
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/future/solving_strategies/schemes/implicit_scheme/BuildWithThreadLocal
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/future/solving_strategies/schemes/implicit_scheme/BuildWithLocalAllocation
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/future/solving_strategies/schemes/implicit_scheme/ImplicitScheme/BuildMassMatrix
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/future/solving_strategies/schemes/implicit_scheme/ImplicitScheme/BuildDampingMatrix
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/future/solving_strategies/schemes/implicit_scheme/ImplicitScheme/BuildMasterSlaveConstraints
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/gid_io/WriteNodalResults
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/gid_io/WriteLocalAxesOnNodes
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/gid_io/WriteNodalFlags
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/gid_io/WriteNodalResultsNonHistorical
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/gid_io/WriteLocalAxesOnNodesNonHistorical
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/gid_io/WriteNodeMesh
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/gid_io/WriteSphereMesh
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/gid_io/WriteCircleMesh
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/gid_io/WriteClusterMesh
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/gid_io/WriteMesh
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/gid_io/PrintFlagsOnGaussPoints
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/gid_io/PrintOnGaussPoints
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/cartesian_mesh_generator_modeler/CartesianMeshGeneratorModeler/GenerateMesh
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/clean_up_problematic_triangles_modeler/CleanUpProblematicTrianglesModeler/SetupModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/duplicate_mesh_modeler/DuplicateMeshModeler/GenerateMesh
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/edge_swapping_2d_modeler/EdgeSwapping2DModeler/Remesh
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/voxel_mesh_generator_modeler/VoxelMeshGeneratorModeler/ApplyColoring
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/voxel_mesh_generator_modeler/VoxelMeshGeneratorModeler/GenerateEntities
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/voxel_mesh_generator_modeler/VoxelMeshGeneratorModeler/ApplyOperations
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/GenerateOctree
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/GenerateNodes
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/CalculateDistance2
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/CalculateDistance
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/GenerateOctree
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/GenerateNodes
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/CalculateDistance2
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/CalculateDistance
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/timer_process/TimerProcess/ExecuteFinalize
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver/ResidualBasedBlockBuilderAndSolver/BuildAndSolve
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver/ResidualBasedBlockBuilderAndSolver/BuildAndSolveLinearizedOnPreviousIteration
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver/ResidualBasedBlockBuilderAndSolver/BuildRHSAndSolve
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver/ResidualBasedBlockBuilderAndSolver/BuildRHS
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver/ResidualBasedBlockBuilderAndSolver/ConstructMasterSlaveConstraintsStructure
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver/ResidualBasedBlockBuilderAndSolver/ConstructMatrixStructure
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver_with_lagrange_multiplier/ResidualBasedBlockBuilderAndSolverWithLagrangeMultiplier/ConstructMasterSlaveConstraintsStructure
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver/ResidualBasedEliminationBuilderAndSolver/BuildAndSolve
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver/ResidualBasedEliminationBuilderAndSolver/ConstructMatrixStructure
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_with_constraints/ResidualBasedEliminationBuilderAndSolverWithConstraints/BuildAndSolveWithConstraints
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_with_constraints/ResidualBasedEliminationBuilderAndSolverWithConstraints/BuildRHSWithConstraints
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_with_constraints/ResidualBasedEliminationBuilderAndSolverWithConstraints/ConstructMatrixStructureWithConstraints
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_with_constraints/ResidualBasedEliminationBuilderAndSolverWithConstraints/ConstructRelationMatrixStructure
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/WriteModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/external_libraries/test_delaunator/validatewithtolerance
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/sources/test_memory_pool/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void Timer::Stop(std::string const& rIntervalName)`

# Calls

- [CreateFullLabel](../../../../functions/utilities/timer/CreateFullLabel.md)
- [pop_back](../../../../functions/containers/pointer_vector_set/PointerVectorSet/pop_back.md)
- [GetStartTime](../../../../functions/utilities/timer/TimerData/GetStartTime.md)

# Called by

- [Build](../../../../functions/future/solving_strategies/schemes/implicit_scheme/ImplicitScheme/Build.md)
- [BuildWithSafeAssemble](../../../../functions/future/solving_strategies/schemes/implicit_scheme/ImplicitScheme/BuildWithSafeAssemble.md)
- [BuildWithThreadLocal](../../../../functions/future/solving_strategies/schemes/implicit_scheme/BuildWithThreadLocal.md)
- [BuildWithLocalAllocation](../../../../functions/future/solving_strategies/schemes/implicit_scheme/BuildWithLocalAllocation.md)
- [BuildMassMatrix](../../../../functions/future/solving_strategies/schemes/implicit_scheme/ImplicitScheme/BuildMassMatrix.md)
- [BuildDampingMatrix](../../../../functions/future/solving_strategies/schemes/implicit_scheme/ImplicitScheme/BuildDampingMatrix.md)
- [BuildMasterSlaveConstraints](../../../../functions/future/solving_strategies/schemes/implicit_scheme/ImplicitScheme/BuildMasterSlaveConstraints.md)
- [WriteNodalResults](../../../../functions/includes/gid_io/WriteNodalResults.md)
- [WriteLocalAxesOnNodes](../../../../functions/includes/gid_io/WriteLocalAxesOnNodes.md)
- [WriteNodalFlags](../../../../functions/includes/gid_io/WriteNodalFlags.md)
- [WriteNodalResultsNonHistorical](../../../../functions/includes/gid_io/WriteNodalResultsNonHistorical.md)
- [WriteLocalAxesOnNodesNonHistorical](../../../../functions/includes/gid_io/WriteLocalAxesOnNodesNonHistorical.md)
- [WriteNodeMesh](../../../../functions/includes/gid_io/WriteNodeMesh.md)
- [WriteSphereMesh](../../../../functions/includes/gid_io/WriteSphereMesh.md)
- [WriteCircleMesh](../../../../functions/includes/gid_io/WriteCircleMesh.md)
- [WriteClusterMesh](../../../../functions/includes/gid_io/WriteClusterMesh.md)
- [WriteMesh](../../../../functions/includes/gid_io/WriteMesh.md)
- [PrintFlagsOnGaussPoints](../../../../functions/includes/gid_io/PrintFlagsOnGaussPoints.md)
- [PrintOnGaussPoints](../../../../functions/includes/gid_io/PrintOnGaussPoints.md)
- [GenerateMesh](../../../../functions/modeler/cartesian_mesh_generator_modeler/CartesianMeshGeneratorModeler/GenerateMesh.md)
- [SetupModelPart](../../../../functions/modeler/clean_up_problematic_triangles_modeler/CleanUpProblematicTrianglesModeler/SetupModelPart.md)
- [GenerateMesh](../../../../functions/modeler/duplicate_mesh_modeler/DuplicateMeshModeler/GenerateMesh.md)
- [Remesh](../../../../functions/modeler/edge_swapping_2d_modeler/EdgeSwapping2DModeler/Remesh.md)
- [ApplyColoring](../../../../functions/modeler/voxel_mesh_generator_modeler/VoxelMeshGeneratorModeler/ApplyColoring.md)
- [GenerateEntities](../../../../functions/modeler/voxel_mesh_generator_modeler/VoxelMeshGeneratorModeler/GenerateEntities.md)
- [ApplyOperations](../../../../functions/modeler/voxel_mesh_generator_modeler/VoxelMeshGeneratorModeler/ApplyOperations.md)
- [GenerateOctree](../../../../functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/GenerateOctree.md)
- [GenerateNodes](../../../../functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/GenerateNodes.md)
- [CalculateDistance2](../../../../functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/CalculateDistance2.md)
- [CalculateDistance](../../../../functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/CalculateDistance.md)
- [GenerateOctree](../../../../functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/GenerateOctree.md)
- [GenerateNodes](../../../../functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/GenerateNodes.md)
- [CalculateDistance2](../../../../functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/CalculateDistance2.md)
- [CalculateDistance](../../../../functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/CalculateDistance.md)
- [ExecuteFinalize](../../../../functions/python_scripts/timer_process/TimerProcess/ExecuteFinalize.md)
- [BuildAndSolve](../../../../functions/solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver/ResidualBasedBlockBuilderAndSolver/BuildAndSolve.md)
- [BuildAndSolveLinearizedOnPreviousIteration](../../../../functions/solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver/ResidualBasedBlockBuilderAndSolver/BuildAndSolveLinearizedOnPreviousIteration.md)
- [BuildRHSAndSolve](../../../../functions/solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver/ResidualBasedBlockBuilderAndSolver/BuildRHSAndSolve.md)
- [BuildRHS](../../../../functions/solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver/ResidualBasedBlockBuilderAndSolver/BuildRHS.md)
- [ConstructMasterSlaveConstraintsStructure](../../../../functions/solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver/ResidualBasedBlockBuilderAndSolver/ConstructMasterSlaveConstraintsStructure.md)
- [ConstructMatrixStructure](../../../../functions/solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver/ResidualBasedBlockBuilderAndSolver/ConstructMatrixStructure.md)
- [ConstructMasterSlaveConstraintsStructure](../../../../functions/solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver_with_lagrange_multiplier/ResidualBasedBlockBuilderAndSolverWithLagrangeMultiplier/ConstructMasterSlaveConstraintsStructure.md)
- [BuildAndSolve](../../../../functions/solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver/ResidualBasedEliminationBuilderAndSolver/BuildAndSolve.md)
- [ConstructMatrixStructure](../../../../functions/solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver/ResidualBasedEliminationBuilderAndSolver/ConstructMatrixStructure.md)
- [BuildAndSolveWithConstraints](../../../../functions/solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_with_constraints/ResidualBasedEliminationBuilderAndSolverWithConstraints/BuildAndSolveWithConstraints.md)
- [BuildRHSWithConstraints](../../../../functions/solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_with_constraints/ResidualBasedEliminationBuilderAndSolverWithConstraints/BuildRHSWithConstraints.md)
- [ConstructMatrixStructureWithConstraints](../../../../functions/solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_with_constraints/ResidualBasedEliminationBuilderAndSolverWithConstraints/ConstructMatrixStructureWithConstraints.md)
- [ConstructRelationMatrixStructure](../../../../functions/solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_with_constraints/ResidualBasedEliminationBuilderAndSolverWithConstraints/ConstructRelationMatrixStructure.md)
- [ReadModelPart](../../../../functions/sources/model_part_io/ModelPartIO/ReadModelPart.md)
- [WriteModelPart](../../../../functions/sources/model_part_io/ModelPartIO/WriteModelPart.md)
- [validatewithtolerance](../../../../functions/tests/cpp_tests/external_libraries/test_delaunator/validatewithtolerance.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/sources/test_memory_pool/KRATOS_TEST_CASE_IN_SUITE.md)