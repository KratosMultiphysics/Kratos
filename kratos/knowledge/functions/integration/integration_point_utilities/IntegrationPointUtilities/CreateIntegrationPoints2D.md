---
type: C++ Method
title: CreateIntegrationPoints2D
resource: integration/integration_point_utilities.cpp#L36-L51
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/integration/integration_info/GetQuadratureMethod
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/integration/integration_point_utilities/IntegrationPointUtilities/CreateIntegrationPoints2DGauss
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/integration/integration_info/IntegrationInfo/GetNumberOfIntegrationPointsPerSpan
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/geometries/brep_surface_on_volume/BrepSurfaceOnVolume/CreateIntegrationPoints
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void IntegrationPointUtilities::CreateIntegrationPoints2D( IntegrationPointsArrayType& rIntegrationPoints, const std::vector<double>& rSpansLocalSpaceU, const std::vector<double>& rSpansLocalSpaceV, const IntegrationInfo& rIntegrationInfo)`

# Calls

- [GetQuadratureMethod](../../../../functions/integration/integration_info/GetQuadratureMethod.md)
- [CreateIntegrationPoints2DGauss](../../../../functions/integration/integration_point_utilities/IntegrationPointUtilities/CreateIntegrationPoints2DGauss.md)
- [GetNumberOfIntegrationPointsPerSpan](../../../../functions/integration/integration_info/IntegrationInfo/GetNumberOfIntegrationPointsPerSpan.md)

# Called by

- [CreateIntegrationPoints](../../../../functions/geometries/brep_surface_on_volume/BrepSurfaceOnVolume/CreateIntegrationPoints.md)