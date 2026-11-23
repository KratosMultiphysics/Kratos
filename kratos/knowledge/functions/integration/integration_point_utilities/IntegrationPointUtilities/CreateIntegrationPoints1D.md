---
type: C++ Method
title: CreateIntegrationPoints1D
resource: integration/integration_point_utilities.cpp#L18-L34
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/integration/integration_info/GetQuadratureMethod
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/integration/integration_point_utilities/IntegrationPointUtilities/CreateIntegrationPoints1DGauss
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/integration/integration_info/IntegrationInfo/GetNumberOfIntegrationPointsPerSpan
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/integration/integration_point_utilities/IntegrationPointUtilities/CreateIntegrationPoints1DGrid
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/geometries/brep_curve/BrepCurve/CreateIntegrationPoints
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/brep_curve_on_surface/BrepCurveOnSurface/CreateIntegrationPoints
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/coupling_geometry/CouplingGeometry/CreateIntegrationPoints
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/nurbs_curve_geometry/NurbsCurveGeometry/CreateIntegrationPoints
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/nurbs_curve_on_surface_geometry/NurbsCurveOnSurfaceGeometry/CreateIntegrationPoints
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void IntegrationPointUtilities::CreateIntegrationPoints1D( IntegrationPointsArrayType& rIntegrationPoints, const std::vector<double>& rSpansLocalSpace, const IntegrationInfo& rIntegrationInfo)`

# Calls

- [GetQuadratureMethod](../../../../functions/integration/integration_info/GetQuadratureMethod.md)
- [CreateIntegrationPoints1DGauss](../../../../functions/integration/integration_point_utilities/IntegrationPointUtilities/CreateIntegrationPoints1DGauss.md)
- [GetNumberOfIntegrationPointsPerSpan](../../../../functions/integration/integration_info/IntegrationInfo/GetNumberOfIntegrationPointsPerSpan.md)
- [CreateIntegrationPoints1DGrid](../../../../functions/integration/integration_point_utilities/IntegrationPointUtilities/CreateIntegrationPoints1DGrid.md)

# Called by

- [CreateIntegrationPoints](../../../../functions/geometries/brep_curve/BrepCurve/CreateIntegrationPoints.md)
- [CreateIntegrationPoints](../../../../functions/geometries/brep_curve_on_surface/BrepCurveOnSurface/CreateIntegrationPoints.md)
- [CreateIntegrationPoints](../../../../functions/geometries/coupling_geometry/CouplingGeometry/CreateIntegrationPoints.md)
- [CreateIntegrationPoints](../../../../functions/geometries/nurbs_curve_geometry/NurbsCurveGeometry/CreateIntegrationPoints.md)
- [CreateIntegrationPoints](../../../../functions/geometries/nurbs_curve_on_surface_geometry/NurbsCurveOnSurfaceGeometry/CreateIntegrationPoints.md)