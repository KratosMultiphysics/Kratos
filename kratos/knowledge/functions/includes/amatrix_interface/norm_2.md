---
type: C++ Function
title: norm_2
resource: includes/amatrix_interface.h#L625-L634
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/constraints/slip_constraint/SlipConstraint/ConstructorHelper
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/edge_based_data_structure/EdgeBasedDataStructure/CalculateEdgeDataValues
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/edge_based_data_structure/EdgeBasedDataStructure/CalculateUnitNormal
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/elements/distance_calculation_element_simplex/DistanceCalculationElementSimplex/CalculateLocalSystem
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/elements/levelset_convection_element_simplex/LevelSetConvectionElementSimplex/CalculateLocalSystem
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/coupling_geometry/CouplingGeometry/SpansLocalSpace
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/coupling_geometry/CouplingGeometry/CreateQuadraturePointGeometries
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/geometry/UnitNormal
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/geometry/CalculateDistance
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/nurbs_curve_on_surface_geometry/NurbsCurveOnSurfaceGeometry/ComputeAxisIntersectionSBM
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/nurbs_curve_on_surface_geometry/NurbsCurveOnSurfaceGeometry/DeterminantOfJacobian
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/nurbs_surface_geometry/NurbsSurfaceGeometry/CalculateEstimatedKnotLengthness
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/nurbs_surface_on_volume_geometry/NurbsSurfaceOnVolumeGeometry/SpansLocalSpace
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/oriented_bounding_box/OrientedBoundingBox<TDim>/OrientedBoundingBox
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/oriented_bounding_box/OrientedBoundingBox<2>/OrientedBoundingBox
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/oriented_bounding_box/OrientedBoundingBox<3>/OrientedBoundingBox
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/plane_3d/Plane3D/Plane3D
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/point/Point/Distance
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/prism_interface_3d_6/PrismInterface3D6/DeterminantOfJacobian
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/quadrature_point_curve_on_surface_geometry/QuadraturePointCurveOnSurfaceGeometry/DeterminantOfJacobian
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/quadrilateral_3d_4/ProjectionPointGlobalToLocalSpace
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/linear_solvers/rayleigh_quotient_iteration_eigenvalue_solver/RayleighQuotientIterationEigenvalueSolver/InitializeSystem
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/clean_up_problematic_triangles_modeler/CleanUpProblematicTrianglesModeler/ComputeDistance
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/operation/compute_surrogate_boundary_data/ComputeSurrogateBoundaryData/Execute
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/operation/find_contacts_in_skin_model_part/FindContactsInSkinModelPart/GetClosestContactColorFromNeighbours
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/apply_constant_vectorvalue_process/ApplyConstantVectorValueProcess/ApplyConstantVectorValueProcess
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/apply_ray_casting_process/ApplyRayCastingProcess<TDim>/CalculateCharacteristicLength
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/assign_scalar_input_to_entities_process/AssignScalarInputToEntitiesProcess<TEntity, THistorical>/ComputeExtrapolationWeight
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_discontinuous_distance_to_skin_process/CalculateDiscontinuousDistanceToSkinProcess<TDim>/CheckIfPointIsRepeated
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_discontinuous_distance_to_skin_process/CalculateDiscontinuousDistanceToSkinProcess<TDim>/ComputeIntersectionNormal
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_discontinuous_distance_to_skin_process/CalculateDiscontinuousDistanceToSkinProcess<TDim>/CorrectDistanceOrientation
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_discontinuous_distance_to_skin_process/CalculateDiscontinuousDistanceToSkinProcess<TDim>/ConvertIntersectionPointToEdgeRatio
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/IsIntersectionNodeOnTetEdge
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/IsNewIntersectionNode
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/IsIntersectionOnCorner
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/CalcSignedDistancesToOneIntNode
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/FindIndexNodesOfTriangle2
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/IntersectionTriangleSegment
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/ComputeDiscontinuousInterpolation
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/IsIntersectionNodeOnTetEdge
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/IsNewIntersectionNode
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/IsIntersectionOnCorner
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/CalcSignedDistancesToMoreThanThreeIntNodes
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/GenerateSkinModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/IntersectionTriangleSegment
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/check_same_modelpart_using_skin_distance_process/CheckSameModelPartUsingSkinDistanceProcess<TDim>/Execute
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/compute_nodal_normal_divergence_process/ComputeNodalNormalDivergenceProcess<THistorical>/GetHistoricalNormalVectorField
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/compute_nodal_normal_divergence_process/ComputeNodalNormalDivergenceProcess<THistorical>/GetNonHistoricalNormalVectorField
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/flux_corrected_transport_convection_process/FluxCorrectedTransportConvectionProcess/CalculateSubStepDeltaTime
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/flux_corrected_transport_convection_process/FluxCorrectedTransportConvectionProcess/CalculateResidual
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/flux_corrected_transport_convection_process/FluxCorrectedTransportConvectionProcess/CalculateAntidiffusiveEdgeContributions
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/levelset_convection_process/LevelSetConvectionProcess/Execute
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/levelset_convection_process/LevelSetConvectionProcess/EvaluateNumberOfSubsteps
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/simple_mortar_mapper_process/SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>/ExecuteImplicitMapping
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/assign_vector_by_direction_process/AssignVectorByDirectionProcess/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/assign_vector_by_direction_to_entity_process/AssignVectorByDirectionToEntityProcess/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/spatial_containers/geometrical_objects_bins/GeometricalObjectsBins/PointIsInsideBoundingBoxWithTolerance
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/geometries/test_geometry/VerifyStrainExactness
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/geometries/test_nurbs_curve_on_surface_geometry/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/geometries/test_nurbs_surface/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/geometries/test_triangle_2d_3/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/geometries/test_triangle_3d_3/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/utilities/test_binbased_fast_point_locator/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/utilities/test_math_utils/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/utilities/test_nearest_point_utilities/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_linear_constraints/TestLinearConstraints/test_constraints
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_vector_interface/TestVectorInterface/test_norm
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/coordinate_transformation_utilities/CalculateRotationOperatorPureShapeSensitivities
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/discont_utils/DiscontinuousShapeFunctionsUtilities/CalculateDiscontinuousShapeFunctions
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/enrichment_utilities/CalculateTetrahedraEnrichedShapeFuncions
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/enrichment_utilities/EnrichmentUtilities/CalculateEnrichedShapeFuncions
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/enrichment_utilities_duplicate_dofs/CalculateTetrahedraEnrichedShapeFuncions
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/exact_mortar_segmentation_utility/ExactMortarIntegrationUtility<3, 3, false>/GetExactIntegration
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/exact_mortar_segmentation_utility/ExactMortarIntegrationUtility<3, 4, false>/GetExactIntegration
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/exact_mortar_segmentation_utility/ExactMortarIntegrationUtility<3, 3, false, 4>/GetExactIntegration
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/exact_mortar_segmentation_utility/ExactMortarIntegrationUtility<3, 4, false, 3>/GetExactIntegration
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/exact_mortar_segmentation_utility/ExactMortarIntegrationUtility<3, 3, true>/GetExactIntegration
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/exact_mortar_segmentation_utility/ExactMortarIntegrationUtility<3, 4, true>/GetExactIntegration
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/exact_mortar_segmentation_utility/ExactMortarIntegrationUtility<3, 3, true, 4>/GetExactIntegration
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/exact_mortar_segmentation_utility/ExactMortarIntegrationUtility<3, 4, true, 3>/GetExactIntegration
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/exact_mortar_segmentation_utility/ExactMortarIntegrationUtility<TDim, TNumNodes, TBelong, TNumNodesMaster>/ComputeAnglesIndexes
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/exact_mortar_segmentation_utility/AnglePoints
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/geometrical_projection_utilities/GeometricalProjectionUtilities/FastMinimalDistanceOnLineWithRadius
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/geometrical_projection_utilities/FastProjectDirection
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/geometrical_projection_utilities/FastProjectOnLine
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/geometrical_projection_utilities/FastProjectOnLine2D
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/geometrical_projection_utilities/ProjectIterativeLine2D
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/geometry_tester/GeometryTesterUtility/VerifyStrainExactness
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/geometry_tester/GeometryTesterUtility/VerfiyShapeFunctionsValues
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/geometry_utilities/GeometryUtils/PointDistanceToLineSegment3D
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/geometry_utilities/CalculateTetrahedraDistances
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/geometry_utilities/CalculateTriangleDistances
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/geometry_utilities/CalculateExactDistancesToPlane
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/intersection_utilities/IntersectionUtilities/ComputeTwoPlaneIntersection
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/intersection_utilities/IntersectionUtilities/ComputeThreePlaneIntersection
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/intersection_utilities/ComputeTriangleLineIntersectionInTheSamePlane
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/intersection_utilities/ComputeTetrahedraLineIntersection
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/math_utils/UnitCrossProduct
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/math_utils/OrthonormalBasisHughesMoeller
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/math_utils/OrthonormalBasisFrisvad
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/math_utils/OrthonormalBasisNaive
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/math_utils/VectorsAngle
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/mls_shape_functions_utility/MLSShapeFunctionsUtility/CalculateKernelDerivative
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/mortar_utilities/GaussPointUnitNormal
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/mortar_utilities/ComputeTangentNodeWithLMAndSlip
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/mortar_utilities/ComputeTangentNodeWithSlip
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/normal_calculation_utils/NormalCalculationUtils/ComputeUnitNormalsFromAreaNormals
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/normal_calculation_utils/CalculateOnSimplex
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/normal_calculation_utils/CalculateOnSimplexLowMemory
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/nurbs_utilities/projection_nurbs_geometry_utilities/NewtonRaphsonCurve
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/nurbs_utilities/projection_nurbs_geometry_utilities/LevenbergMarquardtCurve
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/nurbs_utilities/projection_nurbs_geometry_utilities/NewtonRaphsonSurface
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/nurbs_utilities/projection_nurbs_geometry_utilities/LevenbergMarquardtSurface
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/rbf_shape_functions_utility/RBFShapeFunctionsUtility/CalculateShapeFunctions
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/rbf_shape_functions_utility/RBFShapeFunctionsUtility/CalculateShapeFunctionsAndInterpolation
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/rbf_shape_functions_utility/RBFShapeFunctionsUtility/CalculateInverseMultiquadricShapeParameter
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/rve_periodicity_utility/RVEPeriodicityUtility/AuxiliaryAssignPeriodicity
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
  - target: functions/utilities/signed_distance_calculation_utils/SignedDistanceCalculationUtils/CalculateDistances2D
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/split_tetrahedra_utilities/CalculateSplitTetrahedraShapeFuncions
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/tessellation_utilities/curve_tessellation/CurveTessellation/GetClosestPoint
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`typename TExpressionType::data_type norm_2( AMatrix::MatrixExpression<TExpressionType, TCategory> const& TheExpression)`

# Called by

- [ConstructorHelper](../../../functions/constraints/slip_constraint/SlipConstraint/ConstructorHelper.md)
- [CalculateEdgeDataValues](../../../functions/containers/edge_based_data_structure/EdgeBasedDataStructure/CalculateEdgeDataValues.md)
- [CalculateUnitNormal](../../../functions/containers/edge_based_data_structure/EdgeBasedDataStructure/CalculateUnitNormal.md)
- [CalculateLocalSystem](../../../functions/elements/distance_calculation_element_simplex/DistanceCalculationElementSimplex/CalculateLocalSystem.md)
- [CalculateLocalSystem](../../../functions/elements/levelset_convection_element_simplex/LevelSetConvectionElementSimplex/CalculateLocalSystem.md)
- [SpansLocalSpace](../../../functions/geometries/coupling_geometry/CouplingGeometry/SpansLocalSpace.md)
- [CreateQuadraturePointGeometries](../../../functions/geometries/coupling_geometry/CouplingGeometry/CreateQuadraturePointGeometries.md)
- [UnitNormal](../../../functions/geometries/geometry/UnitNormal.md)
- [CalculateDistance](../../../functions/geometries/geometry/CalculateDistance.md)
- [ComputeAxisIntersectionSBM](../../../functions/geometries/nurbs_curve_on_surface_geometry/NurbsCurveOnSurfaceGeometry/ComputeAxisIntersectionSBM.md)
- [DeterminantOfJacobian](../../../functions/geometries/nurbs_curve_on_surface_geometry/NurbsCurveOnSurfaceGeometry/DeterminantOfJacobian.md)
- [CalculateEstimatedKnotLengthness](../../../functions/geometries/nurbs_surface_geometry/NurbsSurfaceGeometry/CalculateEstimatedKnotLengthness.md)
- [SpansLocalSpace](../../../functions/geometries/nurbs_surface_on_volume_geometry/NurbsSurfaceOnVolumeGeometry/SpansLocalSpace.md)
- [OrientedBoundingBox](../../../functions/geometries/oriented_bounding_box/OrientedBoundingBox<TDim>/OrientedBoundingBox.md)
- [OrientedBoundingBox](../../../functions/geometries/oriented_bounding_box/OrientedBoundingBox<2>/OrientedBoundingBox.md)
- [OrientedBoundingBox](../../../functions/geometries/oriented_bounding_box/OrientedBoundingBox<3>/OrientedBoundingBox.md)
- [Plane3D](../../../functions/geometries/plane_3d/Plane3D/Plane3D.md)
- [Distance](../../../functions/geometries/point/Point/Distance.md)
- [DeterminantOfJacobian](../../../functions/geometries/prism_interface_3d_6/PrismInterface3D6/DeterminantOfJacobian.md)
- [DeterminantOfJacobian](../../../functions/geometries/quadrature_point_curve_on_surface_geometry/QuadraturePointCurveOnSurfaceGeometry/DeterminantOfJacobian.md)
- [ProjectionPointGlobalToLocalSpace](../../../functions/geometries/quadrilateral_3d_4/ProjectionPointGlobalToLocalSpace.md)
- [InitializeSystem](../../../functions/linear_solvers/rayleigh_quotient_iteration_eigenvalue_solver/RayleighQuotientIterationEigenvalueSolver/InitializeSystem.md)
- [ComputeDistance](../../../functions/modeler/clean_up_problematic_triangles_modeler/CleanUpProblematicTrianglesModeler/ComputeDistance.md)
- [Execute](../../../functions/modeler/operation/compute_surrogate_boundary_data/ComputeSurrogateBoundaryData/Execute.md)
- [GetClosestContactColorFromNeighbours](../../../functions/modeler/operation/find_contacts_in_skin_model_part/FindContactsInSkinModelPart/GetClosestContactColorFromNeighbours.md)
- [ApplyConstantVectorValueProcess](../../../functions/processes/apply_constant_vectorvalue_process/ApplyConstantVectorValueProcess/ApplyConstantVectorValueProcess.md)
- [CalculateCharacteristicLength](../../../functions/processes/apply_ray_casting_process/ApplyRayCastingProcess<TDim>/CalculateCharacteristicLength.md)
- [ComputeExtrapolationWeight](../../../functions/processes/assign_scalar_input_to_entities_process/AssignScalarInputToEntitiesProcess<TEntity, THistorical>/ComputeExtrapolationWeight.md)
- [CheckIfPointIsRepeated](../../../functions/processes/calculate_discontinuous_distance_to_skin_process/CalculateDiscontinuousDistanceToSkinProcess<TDim>/CheckIfPointIsRepeated.md)
- [ComputeIntersectionNormal](../../../functions/processes/calculate_discontinuous_distance_to_skin_process/CalculateDiscontinuousDistanceToSkinProcess<TDim>/ComputeIntersectionNormal.md)
- [CorrectDistanceOrientation](../../../functions/processes/calculate_discontinuous_distance_to_skin_process/CalculateDiscontinuousDistanceToSkinProcess<TDim>/CorrectDistanceOrientation.md)
- [ConvertIntersectionPointToEdgeRatio](../../../functions/processes/calculate_discontinuous_distance_to_skin_process/CalculateDiscontinuousDistanceToSkinProcess<TDim>/ConvertIntersectionPointToEdgeRatio.md)
- [IsIntersectionNodeOnTetEdge](../../../functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/IsIntersectionNodeOnTetEdge.md)
- [IsNewIntersectionNode](../../../functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/IsNewIntersectionNode.md)
- [IsIntersectionOnCorner](../../../functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/IsIntersectionOnCorner.md)
- [CalcSignedDistancesToOneIntNode](../../../functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/CalcSignedDistancesToOneIntNode.md)
- [FindIndexNodesOfTriangle2](../../../functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/FindIndexNodesOfTriangle2.md)
- [IntersectionTriangleSegment](../../../functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/IntersectionTriangleSegment.md)
- [ComputeDiscontinuousInterpolation](../../../functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/ComputeDiscontinuousInterpolation.md)
- [IsIntersectionNodeOnTetEdge](../../../functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/IsIntersectionNodeOnTetEdge.md)
- [IsNewIntersectionNode](../../../functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/IsNewIntersectionNode.md)
- [IsIntersectionOnCorner](../../../functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/IsIntersectionOnCorner.md)
- [CalcSignedDistancesToMoreThanThreeIntNodes](../../../functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/CalcSignedDistancesToMoreThanThreeIntNodes.md)
- [GenerateSkinModelPart](../../../functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/GenerateSkinModelPart.md)
- [IntersectionTriangleSegment](../../../functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/IntersectionTriangleSegment.md)
- [Execute](../../../functions/processes/check_same_modelpart_using_skin_distance_process/CheckSameModelPartUsingSkinDistanceProcess<TDim>/Execute.md)
- [GetHistoricalNormalVectorField](../../../functions/processes/compute_nodal_normal_divergence_process/ComputeNodalNormalDivergenceProcess<THistorical>/GetHistoricalNormalVectorField.md)
- [GetNonHistoricalNormalVectorField](../../../functions/processes/compute_nodal_normal_divergence_process/ComputeNodalNormalDivergenceProcess<THistorical>/GetNonHistoricalNormalVectorField.md)
- [CalculateSubStepDeltaTime](../../../functions/processes/flux_corrected_transport_convection_process/FluxCorrectedTransportConvectionProcess/CalculateSubStepDeltaTime.md)
- [CalculateResidual](../../../functions/processes/flux_corrected_transport_convection_process/FluxCorrectedTransportConvectionProcess/CalculateResidual.md)
- [CalculateAntidiffusiveEdgeContributions](../../../functions/processes/flux_corrected_transport_convection_process/FluxCorrectedTransportConvectionProcess/CalculateAntidiffusiveEdgeContributions.md)
- [Execute](../../../functions/processes/levelset_convection_process/LevelSetConvectionProcess/Execute.md)
- [EvaluateNumberOfSubsteps](../../../functions/processes/levelset_convection_process/LevelSetConvectionProcess/EvaluateNumberOfSubsteps.md)
- [ExecuteImplicitMapping](../../../functions/processes/simple_mortar_mapper_process/SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>/ExecuteImplicitMapping.md)
- [__init__](../../../functions/python_scripts/assign_vector_by_direction_process/AssignVectorByDirectionProcess/__init__.md)
- [__init__](../../../functions/python_scripts/assign_vector_by_direction_to_entity_process/AssignVectorByDirectionToEntityProcess/__init__.md)
- [PointIsInsideBoundingBoxWithTolerance](../../../functions/spatial_containers/geometrical_objects_bins/GeometricalObjectsBins/PointIsInsideBoundingBoxWithTolerance.md)
- [VerifyStrainExactness](../../../functions/tests/cpp_tests/geometries/test_geometry/VerifyStrainExactness.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../functions/tests/cpp_tests/geometries/test_nurbs_curve_on_surface_geometry/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../functions/tests/cpp_tests/geometries/test_nurbs_surface/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../functions/tests/cpp_tests/geometries/test_triangle_2d_3/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../functions/tests/cpp_tests/geometries/test_triangle_3d_3/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../functions/tests/cpp_tests/utilities/test_binbased_fast_point_locator/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../functions/tests/cpp_tests/utilities/test_math_utils/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../functions/tests/cpp_tests/utilities/test_nearest_point_utilities/KRATOS_TEST_CASE_IN_SUITE.md)
- [test_constraints](../../../functions/tests/test_linear_constraints/TestLinearConstraints/test_constraints.md)
- [test_norm](../../../functions/tests/test_vector_interface/TestVectorInterface/test_norm.md)
- [CalculateRotationOperatorPureShapeSensitivities](../../../functions/utilities/coordinate_transformation_utilities/CalculateRotationOperatorPureShapeSensitivities.md)
- [CalculateDiscontinuousShapeFunctions](../../../functions/utilities/discont_utils/DiscontinuousShapeFunctionsUtilities/CalculateDiscontinuousShapeFunctions.md)
- [CalculateTetrahedraEnrichedShapeFuncions](../../../functions/utilities/enrichment_utilities/CalculateTetrahedraEnrichedShapeFuncions.md)
- [CalculateEnrichedShapeFuncions](../../../functions/utilities/enrichment_utilities/EnrichmentUtilities/CalculateEnrichedShapeFuncions.md)
- [CalculateTetrahedraEnrichedShapeFuncions](../../../functions/utilities/enrichment_utilities_duplicate_dofs/CalculateTetrahedraEnrichedShapeFuncions.md)
- [GetExactIntegration](../../../functions/utilities/exact_mortar_segmentation_utility/ExactMortarIntegrationUtility<3, 3, false>/GetExactIntegration.md)
- [GetExactIntegration](../../../functions/utilities/exact_mortar_segmentation_utility/ExactMortarIntegrationUtility<3, 4, false>/GetExactIntegration.md)
- [GetExactIntegration](../../../functions/utilities/exact_mortar_segmentation_utility/ExactMortarIntegrationUtility<3, 3, false, 4>/GetExactIntegration.md)
- [GetExactIntegration](../../../functions/utilities/exact_mortar_segmentation_utility/ExactMortarIntegrationUtility<3, 4, false, 3>/GetExactIntegration.md)
- [GetExactIntegration](../../../functions/utilities/exact_mortar_segmentation_utility/ExactMortarIntegrationUtility<3, 3, true>/GetExactIntegration.md)
- [GetExactIntegration](../../../functions/utilities/exact_mortar_segmentation_utility/ExactMortarIntegrationUtility<3, 4, true>/GetExactIntegration.md)
- [GetExactIntegration](../../../functions/utilities/exact_mortar_segmentation_utility/ExactMortarIntegrationUtility<3, 3, true, 4>/GetExactIntegration.md)
- [GetExactIntegration](../../../functions/utilities/exact_mortar_segmentation_utility/ExactMortarIntegrationUtility<3, 4, true, 3>/GetExactIntegration.md)
- [ComputeAnglesIndexes](../../../functions/utilities/exact_mortar_segmentation_utility/ExactMortarIntegrationUtility<TDim, TNumNodes, TBelong, TNumNodesMaster>/ComputeAnglesIndexes.md)
- [AnglePoints](../../../functions/utilities/exact_mortar_segmentation_utility/AnglePoints.md)
- [FastMinimalDistanceOnLineWithRadius](../../../functions/utilities/geometrical_projection_utilities/GeometricalProjectionUtilities/FastMinimalDistanceOnLineWithRadius.md)
- [FastProjectDirection](../../../functions/utilities/geometrical_projection_utilities/FastProjectDirection.md)
- [FastProjectOnLine](../../../functions/utilities/geometrical_projection_utilities/FastProjectOnLine.md)
- [FastProjectOnLine2D](../../../functions/utilities/geometrical_projection_utilities/FastProjectOnLine2D.md)
- [ProjectIterativeLine2D](../../../functions/utilities/geometrical_projection_utilities/ProjectIterativeLine2D.md)
- [VerifyStrainExactness](../../../functions/utilities/geometry_tester/GeometryTesterUtility/VerifyStrainExactness.md)
- [VerfiyShapeFunctionsValues](../../../functions/utilities/geometry_tester/GeometryTesterUtility/VerfiyShapeFunctionsValues.md)
- [PointDistanceToLineSegment3D](../../../functions/utilities/geometry_utilities/GeometryUtils/PointDistanceToLineSegment3D.md)
- [CalculateTetrahedraDistances](../../../functions/utilities/geometry_utilities/CalculateTetrahedraDistances.md)
- [CalculateTriangleDistances](../../../functions/utilities/geometry_utilities/CalculateTriangleDistances.md)
- [CalculateExactDistancesToPlane](../../../functions/utilities/geometry_utilities/CalculateExactDistancesToPlane.md)
- [ComputeTwoPlaneIntersection](../../../functions/utilities/intersection_utilities/IntersectionUtilities/ComputeTwoPlaneIntersection.md)
- [ComputeThreePlaneIntersection](../../../functions/utilities/intersection_utilities/IntersectionUtilities/ComputeThreePlaneIntersection.md)
- [ComputeTriangleLineIntersectionInTheSamePlane](../../../functions/utilities/intersection_utilities/ComputeTriangleLineIntersectionInTheSamePlane.md)
- [ComputeTetrahedraLineIntersection](../../../functions/utilities/intersection_utilities/ComputeTetrahedraLineIntersection.md)
- [UnitCrossProduct](../../../functions/utilities/math_utils/UnitCrossProduct.md)
- [OrthonormalBasisHughesMoeller](../../../functions/utilities/math_utils/OrthonormalBasisHughesMoeller.md)
- [OrthonormalBasisFrisvad](../../../functions/utilities/math_utils/OrthonormalBasisFrisvad.md)
- [OrthonormalBasisNaive](../../../functions/utilities/math_utils/OrthonormalBasisNaive.md)
- [VectorsAngle](../../../functions/utilities/math_utils/VectorsAngle.md)
- [CalculateKernelDerivative](../../../functions/utilities/mls_shape_functions_utility/MLSShapeFunctionsUtility/CalculateKernelDerivative.md)
- [GaussPointUnitNormal](../../../functions/utilities/mortar_utilities/GaussPointUnitNormal.md)
- [ComputeTangentNodeWithLMAndSlip](../../../functions/utilities/mortar_utilities/ComputeTangentNodeWithLMAndSlip.md)
- [ComputeTangentNodeWithSlip](../../../functions/utilities/mortar_utilities/ComputeTangentNodeWithSlip.md)
- [ComputeUnitNormalsFromAreaNormals](../../../functions/utilities/normal_calculation_utils/NormalCalculationUtils/ComputeUnitNormalsFromAreaNormals.md)
- [CalculateOnSimplex](../../../functions/utilities/normal_calculation_utils/CalculateOnSimplex.md)
- [CalculateOnSimplexLowMemory](../../../functions/utilities/normal_calculation_utils/CalculateOnSimplexLowMemory.md)
- [NewtonRaphsonCurve](../../../functions/utilities/nurbs_utilities/projection_nurbs_geometry_utilities/NewtonRaphsonCurve.md)
- [LevenbergMarquardtCurve](../../../functions/utilities/nurbs_utilities/projection_nurbs_geometry_utilities/LevenbergMarquardtCurve.md)
- [NewtonRaphsonSurface](../../../functions/utilities/nurbs_utilities/projection_nurbs_geometry_utilities/NewtonRaphsonSurface.md)
- [LevenbergMarquardtSurface](../../../functions/utilities/nurbs_utilities/projection_nurbs_geometry_utilities/LevenbergMarquardtSurface.md)
- [CalculateShapeFunctions](../../../functions/utilities/rbf_shape_functions_utility/RBFShapeFunctionsUtility/CalculateShapeFunctions.md)
- [CalculateShapeFunctionsAndInterpolation](../../../functions/utilities/rbf_shape_functions_utility/RBFShapeFunctionsUtility/CalculateShapeFunctionsAndInterpolation.md)
- [CalculateInverseMultiquadricShapeParameter](../../../functions/utilities/rbf_shape_functions_utility/RBFShapeFunctionsUtility/CalculateInverseMultiquadricShapeParameter.md)
- [AuxiliaryAssignPeriodicity](../../../functions/utilities/rve_periodicity_utility/RVEPeriodicityUtility/AuxiliaryAssignPeriodicity.md)
- [CalculateGradientBasedConformingExtensionBasis](../../../functions/utilities/shifted_boundary_meshless_interface_utility/ShiftedBoundaryMeshlessInterfaceUtility/CalculateGradientBasedConformingExtensionBasis.md)
- [CalculateMeshlessBasedConformingExtensionBasis](../../../functions/utilities/shifted_boundary_meshless_interface_utility/ShiftedBoundaryMeshlessInterfaceUtility/CalculateMeshlessBasedConformingExtensionBasis.md)
- [CalculateMeshlessBasedNonConformingExtensionBasis](../../../functions/utilities/shifted_boundary_meshless_interface_utility/ShiftedBoundaryMeshlessInterfaceUtility/CalculateMeshlessBasedNonConformingExtensionBasis.md)
- [CalculateDistances2D](../../../functions/utilities/signed_distance_calculation_utils/SignedDistanceCalculationUtils/CalculateDistances2D.md)
- [CalculateSplitTetrahedraShapeFuncions](../../../functions/utilities/split_tetrahedra_utilities/CalculateSplitTetrahedraShapeFuncions.md)
- [GetClosestPoint](../../../functions/utilities/tessellation_utilities/curve_tessellation/CurveTessellation/GetClosestPoint.md)