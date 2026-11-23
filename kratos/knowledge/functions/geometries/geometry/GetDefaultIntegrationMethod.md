---
type: C++ Function
title: GetDefaultIntegrationMethod
resource: geometries/geometry.h#L2066-L2069
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/conditions/mesh_condition/MeshCondition/CalculateOnIntegrationPoints
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/elements/mesh_element/MeshElement/CalculateOnIntegrationPoints
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/geometry/GetDefaultIntegrationInfo
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/geometry/GlobalCoordinates
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/nurbs_curve_geometry/NurbsCurveGeometry/CreateQuadraturePointGeometries
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/nurbs_curve_on_surface_geometry/NurbsCurveOnSurfaceGeometry/CreateQuadraturePointGeometries
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/nurbs_curve_on_surface_geometry/NurbsCurveOnSurfaceGeometry/CreateQuadraturePointGeometriesSBM
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/nurbs_surface_geometry/NurbsSurfaceGeometry/CreateQuadraturePointGeometries
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/nurbs_surface_on_volume_geometry/NurbsSurfaceOnVolumeGeometry/CreateQuadraturePointGeometries
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/nurbs_surface_on_volume_geometry/NurbsSurfaceOnVolumeGeometry/CreateQuadraturePointGeometriesSbm
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/nurbs_volume_geometry/NurbsVolumeGeometry/CreateQuadraturePointGeometries
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/surface_in_nurbs_volume_geometry/SurfaceInNurbsVolumeGeometry/CreateIntegrationPoints
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/surface_in_nurbs_volume_geometry/SurfaceInNurbsVolumeGeometry/CreateQuadraturePointGeometries
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/condition/Condition/GetIntegrationMethod
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/element/Element/GetIntegrationMethod
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_nodal_area_process/CalculateNodalAreaProcess<THistorical>/Execute
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/compute_nodal_gradient_process/ComputeNodalGradientProcess<TOutputHistorical>/ComputeElementalContributionsAndVolume
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/compute_nodal_normal_divergence_process/ComputeNodalNormalDivergenceProcess<THistorical>/Execute
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python/add_geometries_to_python/AddGeometriesToPython
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tensor_adaptors/geometries_tensor_adaptor/GeometriesTensorAdaptor/GetIntegrationMethod
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/geometries/test_non_rectangular_jacobian/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/geometries/test_nurbs_curve_on_surface_geometry/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/geometries/test_surface_in_nurbs_volume/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/sources/test_serializer_geometry_data/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/sources/test_shape_function_derivatives/TestAllShapeFunctionsLocalGradients
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/utilities/test_geometry_utils/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_geometry_data/TestGeometryData/test_geometry_data
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_utilities/test_element/TestElement/Initialize
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/integration_utilities/GetIntegrationMethodForExactMassMatrixEvaluation
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/integration_utilities/ComputeArea2DGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/integration_utilities/ComputeArea2DGeometryDerivative
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/integration_utilities/ComputeVolume3DGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/integration_utilities/ComputeVolume3DGeometryDerivative
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/quadrature_points_utility/CreateQuadraturePointsUtility/CreateFromLocalCoordinates
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/quadrature_points_utility/CreateQuadraturePointsUtility/Create
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/quadrature_points_utility/CreateQuadraturePointsUtility/UpdateFromLocalCoordinates
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`IntegrationMethod GetDefaultIntegrationMethod() const`

# Called by

- [CalculateOnIntegrationPoints](../../../functions/conditions/mesh_condition/MeshCondition/CalculateOnIntegrationPoints.md)
- [CalculateOnIntegrationPoints](../../../functions/elements/mesh_element/MeshElement/CalculateOnIntegrationPoints.md)
- [GetDefaultIntegrationInfo](../../../functions/geometries/geometry/GetDefaultIntegrationInfo.md)
- [GlobalCoordinates](../../../functions/geometries/geometry/GlobalCoordinates.md)
- [CreateQuadraturePointGeometries](../../../functions/geometries/nurbs_curve_geometry/NurbsCurveGeometry/CreateQuadraturePointGeometries.md)
- [CreateQuadraturePointGeometries](../../../functions/geometries/nurbs_curve_on_surface_geometry/NurbsCurveOnSurfaceGeometry/CreateQuadraturePointGeometries.md)
- [CreateQuadraturePointGeometriesSBM](../../../functions/geometries/nurbs_curve_on_surface_geometry/NurbsCurveOnSurfaceGeometry/CreateQuadraturePointGeometriesSBM.md)
- [CreateQuadraturePointGeometries](../../../functions/geometries/nurbs_surface_geometry/NurbsSurfaceGeometry/CreateQuadraturePointGeometries.md)
- [CreateQuadraturePointGeometries](../../../functions/geometries/nurbs_surface_on_volume_geometry/NurbsSurfaceOnVolumeGeometry/CreateQuadraturePointGeometries.md)
- [CreateQuadraturePointGeometriesSbm](../../../functions/geometries/nurbs_surface_on_volume_geometry/NurbsSurfaceOnVolumeGeometry/CreateQuadraturePointGeometriesSbm.md)
- [CreateQuadraturePointGeometries](../../../functions/geometries/nurbs_volume_geometry/NurbsVolumeGeometry/CreateQuadraturePointGeometries.md)
- [CreateIntegrationPoints](../../../functions/geometries/surface_in_nurbs_volume_geometry/SurfaceInNurbsVolumeGeometry/CreateIntegrationPoints.md)
- [CreateQuadraturePointGeometries](../../../functions/geometries/surface_in_nurbs_volume_geometry/SurfaceInNurbsVolumeGeometry/CreateQuadraturePointGeometries.md)
- [GetIntegrationMethod](../../../functions/includes/condition/Condition/GetIntegrationMethod.md)
- [GetIntegrationMethod](../../../functions/includes/element/Element/GetIntegrationMethod.md)
- [Execute](../../../functions/processes/calculate_nodal_area_process/CalculateNodalAreaProcess<THistorical>/Execute.md)
- [ComputeElementalContributionsAndVolume](../../../functions/processes/compute_nodal_gradient_process/ComputeNodalGradientProcess<TOutputHistorical>/ComputeElementalContributionsAndVolume.md)
- [Execute](../../../functions/processes/compute_nodal_normal_divergence_process/ComputeNodalNormalDivergenceProcess<THistorical>/Execute.md)
- [AddGeometriesToPython](../../../functions/python/add_geometries_to_python/AddGeometriesToPython.md)
- [GetIntegrationMethod](../../../functions/tensor_adaptors/geometries_tensor_adaptor/GeometriesTensorAdaptor/GetIntegrationMethod.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../functions/tests/cpp_tests/geometries/test_non_rectangular_jacobian/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../functions/tests/cpp_tests/geometries/test_nurbs_curve_on_surface_geometry/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../functions/tests/cpp_tests/geometries/test_surface_in_nurbs_volume/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../functions/tests/cpp_tests/sources/test_serializer_geometry_data/KRATOS_TEST_CASE_IN_SUITE.md)
- [TestAllShapeFunctionsLocalGradients](../../../functions/tests/cpp_tests/sources/test_shape_function_derivatives/TestAllShapeFunctionsLocalGradients.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../functions/tests/cpp_tests/utilities/test_geometry_utils/KRATOS_TEST_CASE_IN_SUITE.md)
- [test_geometry_data](../../../functions/tests/test_geometry_data/TestGeometryData/test_geometry_data.md)
- [Initialize](../../../functions/tests/test_utilities/test_element/TestElement/Initialize.md)
- [GetIntegrationMethodForExactMassMatrixEvaluation](../../../functions/utilities/integration_utilities/GetIntegrationMethodForExactMassMatrixEvaluation.md)
- [ComputeArea2DGeometry](../../../functions/utilities/integration_utilities/ComputeArea2DGeometry.md)
- [ComputeArea2DGeometryDerivative](../../../functions/utilities/integration_utilities/ComputeArea2DGeometryDerivative.md)
- [ComputeVolume3DGeometry](../../../functions/utilities/integration_utilities/ComputeVolume3DGeometry.md)
- [ComputeVolume3DGeometryDerivative](../../../functions/utilities/integration_utilities/ComputeVolume3DGeometryDerivative.md)
- [CreateFromLocalCoordinates](../../../functions/utilities/quadrature_points_utility/CreateQuadraturePointsUtility/CreateFromLocalCoordinates.md)
- [Create](../../../functions/utilities/quadrature_points_utility/CreateQuadraturePointsUtility/Create.md)
- [UpdateFromLocalCoordinates](../../../functions/utilities/quadrature_points_utility/CreateQuadraturePointsUtility/UpdateFromLocalCoordinates.md)