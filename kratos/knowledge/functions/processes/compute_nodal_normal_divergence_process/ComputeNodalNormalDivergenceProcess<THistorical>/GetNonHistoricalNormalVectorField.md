---
type: C++ Method
title: GetNonHistoricalNormalVectorField
resource: processes/compute_nodal_normal_divergence_process.cpp#L308-L318
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/amatrix_interface/norm_2
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`array_1d<double,3> ComputeNodalNormalDivergenceProcess<THistorical>::GetNonHistoricalNormalVectorField( const Node& rNode, const Variable<array_1d<double,3>>& rVariable)`

# Calls

- [norm_2](../../../../functions/includes/amatrix_interface/norm_2.md)