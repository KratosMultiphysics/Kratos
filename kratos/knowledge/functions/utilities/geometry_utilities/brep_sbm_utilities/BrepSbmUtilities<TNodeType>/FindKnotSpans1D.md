---
type: C++ Method
title: FindKnotSpans1D
resource: utilities/geometry_utilities/brep_sbm_utilities.cpp#L274-L286
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/utilities/geometry_utilities/brep_sbm_utilities/BrepSbmUtilities<TNodeType>/CreateBrepSurfaceSbmIntegrationPoints
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/geometry_utilities/brep_sbm_utilities/BrepSbmUtilities<TNodeType>/CreateBrepVolumeSbmIntegrationPoints
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`int BrepSbmUtilities<TNodeType>::FindKnotSpans1D( const std::vector<double>& rSpans, const double coord)`

# Called by

- [CreateBrepSurfaceSbmIntegrationPoints](../../../../../functions/utilities/geometry_utilities/brep_sbm_utilities/BrepSbmUtilities<TNodeType>/CreateBrepSurfaceSbmIntegrationPoints.md)
- [CreateBrepVolumeSbmIntegrationPoints](../../../../../functions/utilities/geometry_utilities/brep_sbm_utilities/BrepSbmUtilities<TNodeType>/CreateBrepVolumeSbmIntegrationPoints.md)