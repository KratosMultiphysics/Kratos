---
type: C++ Function
title: prod
resource: includes/amatrix_interface.h#L568-L578
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/elements/distance_calculation_element_simplex/DistanceCalculationElementSimplex/CalculateLocalSystem
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/elements/levelset_convection_element_simplex/LevelSetConvectionElementSimplex/CalculateLocalSystem
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/elements/levelset_convection_element_simplex_algebraic_stabilization/LevelSetConvectionElementSimplexAlgebraicStabilization<TDim, TNumNodes>/CalculateLocalSystem
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/geometry/ShapeFunctionsIntegrationPointsGradients
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/nurbs_surface_on_volume_geometry/NurbsSurfaceOnVolumeGeometry/GlobalSpaceDerivatives
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/oriented_bounding_box/OrientedBoundingBox<TDim>/RotateNode3D
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/quadrature_point_surface_in_volume_geometry/QuadraturePointSurfaceInVolumeGeometry/Normal
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/quadrature_point_surface_in_volume_geometry/QuadraturePointSurfaceInVolumeGeometry/DeterminantOfJacobian
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/surface_in_nurbs_volume_geometry/SurfaceInNurbsVolumeGeometry/GlobalSpaceDerivatives
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/constitutive_law/AddInitialDeformationGradientMatrixContribution
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/mortar_classes/MortarOperator/ComputePOperator
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/mortar_classes/DualLagrangeMultiplierOperators/CalculateAe
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/mortar_classes/DualLagrangeMultiplierOperatorsWithDerivatives/CalculateDeltaAe
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modified_shape_functions/modified_shape_functions/ModifiedShapeFunctions/ComputeValuesOnOneSide
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modified_shape_functions/modified_shape_functions/ModifiedShapeFunctions/ComputeFaceValuesOnOneSide
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/ComputeDiscontinuousInterpolation
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/CalcSignedDistancesToMoreThanThreeIntNodes
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/GenerateSkinModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/compute_nodal_gradient_process/ComputeNodalGradientProcess<TOutputHistorical>/ComputeElementalContributionsAndVolume
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/simple_mortar_mapper_process/SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>/AssemblyMortarOperators
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/simple_mortar_mapper_process/SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>/ComputeResidualMatrix
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python/add_matrix_to_python/CreateMatrixInterface
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python/add_matrix_to_python/AddMatrixToPython
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/CalculateEntityResidualLocalContributions
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/CalculateEntityTimeSchemeContributions
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/CalculateEntityAuxiliaryVariableContributions
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/schemes/residual_based_adjoint_static_scheme/ResidualBasedAdjointStaticScheme/CalculateSystemContributions
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/schemes/residual_based_adjoint_static_scheme/ResidualBasedAdjointStaticScheme/CalculateRHSContribution
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/schemes/residual_based_adjoint_steady_scheme/ResidualBasedAdjointSteadyScheme/CalculateSystemContributions
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/schemes/residual_based_bdf_scheme/TemplateAddDynamicsToRHS
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/schemes/residual_based_bossak_displacement_scheme/ResidualBasedBossakDisplacementScheme/AddDynamicsToRHS
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/schemes/residual_based_pseudo_static_displacement_scheme/ResidualBasedPseudoStaticDisplacementScheme/AddDynamicsToRHS
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/schemes/sensitivity_builder_scheme/CalculateLocalSensitivity
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/constitutive_law/ConstitutiveLaw/ContraVariantPushForward
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/constitutive_law/ConstitutiveLaw/ContraVariantPullBack
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/constitutive_law/ConstitutiveLaw/CoVariantPushForward
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/constitutive_law/ConstitutiveLaw/CoVariantPullBack
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/geometries/test_geometry/VerifyStrainExactness
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/strategies/strategies/explicit_strategy_testing_utils/LogFittingSlope
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/utilities/test_dense_householder_qr_decomposition/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/utilities/test_line_geometrical_sensitivity/GetJacobian
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/utilities/test_line_geometrical_sensitivity/IntegrationPointWeight
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/utilities/test_math_utils/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/utilities/test_svd_utils/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_utilities/test_laplacian_element/TestLaplacianElement/CalculateLocalSystem
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/beam_math_utilities/InverseCayleyTransform
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/beam_math_utilities/BeamMathUtils/MapLocalToGlobal2D
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/beam_math_utilities/BeamMathUtils/MapLocalToGlobal3D
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/beam_math_utilities/BeamMathUtils/CheckOrthogonality
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/constraint_utilities/PreComputeExplicitConstraintConstribution
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/constraint_utilities/PreComputeExplicitConstraintMassAndInertia
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/coordinate_transformation_utilities/RotateVelocities
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/coordinate_transformation_utilities/RecoverVelocities
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/coordinate_transformation_utilities/RotateAuxPure
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/coordinate_transformation_utilities/RotateRHSAux
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
  - target: functions/utilities/enrichment_utilities/EnrichmentUtilities/CalculateEnrichedShapeFuncionsExtended
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/enrichment_utilities/EnrichmentUtilities/CalculateEnrichedShapeFuncionsInLocalAxis
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/enrichment_utilities_duplicate_dofs/CalculateTetrahedraEnrichedShapeFuncions
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/function_parser_utility/GenericFunctionUtility/RotateAndCallFunction
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/geometrical_projection_utilities/ProjectIterativeLine2D
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/geometrical_sensitivity_utility/GeometricalSensitivityUtility/CalculateSensitivity
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/geometry_tester/GeometryTesterUtility/VerifyStrainExactness
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/geometry_utilities/GeometryUtils/ShapeFunctionsSecondDerivativesTransformOnIntegrationPoint
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/geometry_utilities/CalculateExactDistancesToPlane
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/geometry_utilities/ShapeFunctionsGradients
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/geometry_utilities/DeformationGradient
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/line_sensitivity_utility/LineSensitivityUtility/CalculateSensitivity
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/line_sensitivity_utility/LineSensitivityUtility/Initialize
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/math_utils/GeneralizedInvertMatrix
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/math_utils/GeneralizedDet
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/math_utils/CalculateExponentialOfMatrix
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/mortar_utilities/RotatePoint
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/qef_utility/QuadraticErrorFunction/UpdateATAATB
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/qef_utility/QuadraticErrorFunction/ComputeQuadraticErrorFunctionPoint
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/rve_periodicity_utility/RVEPeriodicityUtility/AuxiliaryAssignPeriodicity
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/signed_distance_calculation_utils/SignedDistanceCalculationUtils/CalculateDistances2D
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/split_tetrahedra_utilities/CalculateSplitTetrahedraShapeFuncions
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/svd_utils/SVDUtils/SingularValueDecomposition2x2
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/svd_utils/SVDUtils/SingularValueDecomposition2x2Symmetric
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/svd_utils/ComputeJacobiIterations
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`AMatrix::MatrixProductExpression<TExpression1Type, TExpression2Type> prod( AMatrix::MatrixExpression<TExpression1Type, TCategory1> const& First, AMatrix::MatrixExpression<TExpression2Type, TCategory2> const& Second)`

# Called by

- [CalculateLocalSystem](../../../functions/elements/distance_calculation_element_simplex/DistanceCalculationElementSimplex/CalculateLocalSystem.md)
- [CalculateLocalSystem](../../../functions/elements/levelset_convection_element_simplex/LevelSetConvectionElementSimplex/CalculateLocalSystem.md)
- [CalculateLocalSystem](../../../functions/elements/levelset_convection_element_simplex_algebraic_stabilization/LevelSetConvectionElementSimplexAlgebraicStabilization<TDim, TNumNodes>/CalculateLocalSystem.md)
- [ShapeFunctionsIntegrationPointsGradients](../../../functions/geometries/geometry/ShapeFunctionsIntegrationPointsGradients.md)
- [GlobalSpaceDerivatives](../../../functions/geometries/nurbs_surface_on_volume_geometry/NurbsSurfaceOnVolumeGeometry/GlobalSpaceDerivatives.md)
- [RotateNode3D](../../../functions/geometries/oriented_bounding_box/OrientedBoundingBox<TDim>/RotateNode3D.md)
- [Normal](../../../functions/geometries/quadrature_point_surface_in_volume_geometry/QuadraturePointSurfaceInVolumeGeometry/Normal.md)
- [DeterminantOfJacobian](../../../functions/geometries/quadrature_point_surface_in_volume_geometry/QuadraturePointSurfaceInVolumeGeometry/DeterminantOfJacobian.md)
- [GlobalSpaceDerivatives](../../../functions/geometries/surface_in_nurbs_volume_geometry/SurfaceInNurbsVolumeGeometry/GlobalSpaceDerivatives.md)
- [AddInitialDeformationGradientMatrixContribution](../../../functions/includes/constitutive_law/AddInitialDeformationGradientMatrixContribution.md)
- [ComputePOperator](../../../functions/includes/mortar_classes/MortarOperator/ComputePOperator.md)
- [CalculateAe](../../../functions/includes/mortar_classes/DualLagrangeMultiplierOperators/CalculateAe.md)
- [CalculateDeltaAe](../../../functions/includes/mortar_classes/DualLagrangeMultiplierOperatorsWithDerivatives/CalculateDeltaAe.md)
- [ComputeValuesOnOneSide](../../../functions/modified_shape_functions/modified_shape_functions/ModifiedShapeFunctions/ComputeValuesOnOneSide.md)
- [ComputeFaceValuesOnOneSide](../../../functions/modified_shape_functions/modified_shape_functions/ModifiedShapeFunctions/ComputeFaceValuesOnOneSide.md)
- [ComputeDiscontinuousInterpolation](../../../functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/ComputeDiscontinuousInterpolation.md)
- [CalcSignedDistancesToMoreThanThreeIntNodes](../../../functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/CalcSignedDistancesToMoreThanThreeIntNodes.md)
- [GenerateSkinModelPart](../../../functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/GenerateSkinModelPart.md)
- [ComputeElementalContributionsAndVolume](../../../functions/processes/compute_nodal_gradient_process/ComputeNodalGradientProcess<TOutputHistorical>/ComputeElementalContributionsAndVolume.md)
- [AssemblyMortarOperators](../../../functions/processes/simple_mortar_mapper_process/SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>/AssemblyMortarOperators.md)
- [ComputeResidualMatrix](../../../functions/processes/simple_mortar_mapper_process/SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>/ComputeResidualMatrix.md)
- [CreateMatrixInterface](../../../functions/python/add_matrix_to_python/CreateMatrixInterface.md)
- [AddMatrixToPython](../../../functions/python/add_matrix_to_python/AddMatrixToPython.md)
- [CalculateEntityResidualLocalContributions](../../../functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/CalculateEntityResidualLocalContributions.md)
- [CalculateEntityTimeSchemeContributions](../../../functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/CalculateEntityTimeSchemeContributions.md)
- [CalculateEntityAuxiliaryVariableContributions](../../../functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/CalculateEntityAuxiliaryVariableContributions.md)
- [CalculateSystemContributions](../../../functions/solving_strategies/schemes/residual_based_adjoint_static_scheme/ResidualBasedAdjointStaticScheme/CalculateSystemContributions.md)
- [CalculateRHSContribution](../../../functions/solving_strategies/schemes/residual_based_adjoint_static_scheme/ResidualBasedAdjointStaticScheme/CalculateRHSContribution.md)
- [CalculateSystemContributions](../../../functions/solving_strategies/schemes/residual_based_adjoint_steady_scheme/ResidualBasedAdjointSteadyScheme/CalculateSystemContributions.md)
- [TemplateAddDynamicsToRHS](../../../functions/solving_strategies/schemes/residual_based_bdf_scheme/TemplateAddDynamicsToRHS.md)
- [AddDynamicsToRHS](../../../functions/solving_strategies/schemes/residual_based_bossak_displacement_scheme/ResidualBasedBossakDisplacementScheme/AddDynamicsToRHS.md)
- [AddDynamicsToRHS](../../../functions/solving_strategies/schemes/residual_based_pseudo_static_displacement_scheme/ResidualBasedPseudoStaticDisplacementScheme/AddDynamicsToRHS.md)
- [CalculateLocalSensitivity](../../../functions/solving_strategies/schemes/sensitivity_builder_scheme/CalculateLocalSensitivity.md)
- [ContraVariantPushForward](../../../functions/sources/constitutive_law/ConstitutiveLaw/ContraVariantPushForward.md)
- [ContraVariantPullBack](../../../functions/sources/constitutive_law/ConstitutiveLaw/ContraVariantPullBack.md)
- [CoVariantPushForward](../../../functions/sources/constitutive_law/ConstitutiveLaw/CoVariantPushForward.md)
- [CoVariantPullBack](../../../functions/sources/constitutive_law/ConstitutiveLaw/CoVariantPullBack.md)
- [VerifyStrainExactness](../../../functions/tests/cpp_tests/geometries/test_geometry/VerifyStrainExactness.md)
- [LogFittingSlope](../../../functions/tests/cpp_tests/strategies/strategies/explicit_strategy_testing_utils/LogFittingSlope.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../functions/tests/cpp_tests/utilities/test_dense_householder_qr_decomposition/KRATOS_TEST_CASE_IN_SUITE.md)
- [GetJacobian](../../../functions/tests/cpp_tests/utilities/test_line_geometrical_sensitivity/GetJacobian.md)
- [IntegrationPointWeight](../../../functions/tests/cpp_tests/utilities/test_line_geometrical_sensitivity/IntegrationPointWeight.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../functions/tests/cpp_tests/utilities/test_math_utils/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../functions/tests/cpp_tests/utilities/test_svd_utils/KRATOS_TEST_CASE_IN_SUITE.md)
- [CalculateLocalSystem](../../../functions/tests/test_utilities/test_laplacian_element/TestLaplacianElement/CalculateLocalSystem.md)
- [InverseCayleyTransform](../../../functions/utilities/beam_math_utilities/InverseCayleyTransform.md)
- [MapLocalToGlobal2D](../../../functions/utilities/beam_math_utilities/BeamMathUtils/MapLocalToGlobal2D.md)
- [MapLocalToGlobal3D](../../../functions/utilities/beam_math_utilities/BeamMathUtils/MapLocalToGlobal3D.md)
- [CheckOrthogonality](../../../functions/utilities/beam_math_utilities/BeamMathUtils/CheckOrthogonality.md)
- [PreComputeExplicitConstraintConstribution](../../../functions/utilities/constraint_utilities/PreComputeExplicitConstraintConstribution.md)
- [PreComputeExplicitConstraintMassAndInertia](../../../functions/utilities/constraint_utilities/PreComputeExplicitConstraintMassAndInertia.md)
- [RotateVelocities](../../../functions/utilities/coordinate_transformation_utilities/RotateVelocities.md)
- [RecoverVelocities](../../../functions/utilities/coordinate_transformation_utilities/RecoverVelocities.md)
- [RotateAuxPure](../../../functions/utilities/coordinate_transformation_utilities/RotateAuxPure.md)
- [RotateRHSAux](../../../functions/utilities/coordinate_transformation_utilities/RotateRHSAux.md)
- [CalculateDiscontinuousShapeFunctions](../../../functions/utilities/discont_utils/DiscontinuousShapeFunctionsUtilities/CalculateDiscontinuousShapeFunctions.md)
- [CalculateTetrahedraEnrichedShapeFuncions](../../../functions/utilities/enrichment_utilities/CalculateTetrahedraEnrichedShapeFuncions.md)
- [CalculateEnrichedShapeFuncions](../../../functions/utilities/enrichment_utilities/EnrichmentUtilities/CalculateEnrichedShapeFuncions.md)
- [CalculateEnrichedShapeFuncionsExtended](../../../functions/utilities/enrichment_utilities/EnrichmentUtilities/CalculateEnrichedShapeFuncionsExtended.md)
- [CalculateEnrichedShapeFuncionsInLocalAxis](../../../functions/utilities/enrichment_utilities/EnrichmentUtilities/CalculateEnrichedShapeFuncionsInLocalAxis.md)
- [CalculateTetrahedraEnrichedShapeFuncions](../../../functions/utilities/enrichment_utilities_duplicate_dofs/CalculateTetrahedraEnrichedShapeFuncions.md)
- [RotateAndCallFunction](../../../functions/utilities/function_parser_utility/GenericFunctionUtility/RotateAndCallFunction.md)
- [ProjectIterativeLine2D](../../../functions/utilities/geometrical_projection_utilities/ProjectIterativeLine2D.md)
- [CalculateSensitivity](../../../functions/utilities/geometrical_sensitivity_utility/GeometricalSensitivityUtility/CalculateSensitivity.md)
- [VerifyStrainExactness](../../../functions/utilities/geometry_tester/GeometryTesterUtility/VerifyStrainExactness.md)
- [ShapeFunctionsSecondDerivativesTransformOnIntegrationPoint](../../../functions/utilities/geometry_utilities/GeometryUtils/ShapeFunctionsSecondDerivativesTransformOnIntegrationPoint.md)
- [CalculateExactDistancesToPlane](../../../functions/utilities/geometry_utilities/CalculateExactDistancesToPlane.md)
- [ShapeFunctionsGradients](../../../functions/utilities/geometry_utilities/ShapeFunctionsGradients.md)
- [DeformationGradient](../../../functions/utilities/geometry_utilities/DeformationGradient.md)
- [CalculateSensitivity](../../../functions/utilities/line_sensitivity_utility/LineSensitivityUtility/CalculateSensitivity.md)
- [Initialize](../../../functions/utilities/line_sensitivity_utility/LineSensitivityUtility/Initialize.md)
- [GeneralizedInvertMatrix](../../../functions/utilities/math_utils/GeneralizedInvertMatrix.md)
- [GeneralizedDet](../../../functions/utilities/math_utils/GeneralizedDet.md)
- [CalculateExponentialOfMatrix](../../../functions/utilities/math_utils/CalculateExponentialOfMatrix.md)
- [RotatePoint](../../../functions/utilities/mortar_utilities/RotatePoint.md)
- [UpdateATAATB](../../../functions/utilities/qef_utility/QuadraticErrorFunction/UpdateATAATB.md)
- [ComputeQuadraticErrorFunctionPoint](../../../functions/utilities/qef_utility/QuadraticErrorFunction/ComputeQuadraticErrorFunctionPoint.md)
- [AuxiliaryAssignPeriodicity](../../../functions/utilities/rve_periodicity_utility/RVEPeriodicityUtility/AuxiliaryAssignPeriodicity.md)
- [CalculateDistances2D](../../../functions/utilities/signed_distance_calculation_utils/SignedDistanceCalculationUtils/CalculateDistances2D.md)
- [CalculateSplitTetrahedraShapeFuncions](../../../functions/utilities/split_tetrahedra_utilities/CalculateSplitTetrahedraShapeFuncions.md)
- [SingularValueDecomposition2x2](../../../functions/utilities/svd_utils/SVDUtils/SingularValueDecomposition2x2.md)
- [SingularValueDecomposition2x2Symmetric](../../../functions/utilities/svd_utils/SVDUtils/SingularValueDecomposition2x2Symmetric.md)
- [ComputeJacobiIterations](../../../functions/utilities/svd_utils/ComputeJacobiIterations.md)