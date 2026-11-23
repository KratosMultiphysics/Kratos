---
type: C++ Method
title: AppendIdsAndWeights
resource: utilities/rve_periodicity_utility.cpp#L144-L168
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/utilities/rve_periodicity_utility/RVEPeriodicityUtility/Finalize
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void RVEPeriodicityUtility::AppendIdsAndWeights( std::map<IndexType, DataTupletype>& rAux, const IndexType MasterId, const double MasterWeight, std::vector<IndexType>& rFinalMastersIds, std::vector<double>& rFinalMastersWeights, Vector& rFinalT)`

# Called by

- [Finalize](../../../../functions/utilities/rve_periodicity_utility/RVEPeriodicityUtility/Finalize.md)