---
type: C++ Function
title: RecoverVelocities
resource: utilities/coordinate_transformation_utilities.h#L653-L687
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/coordinate_transformation_utilities/IsSlip
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/amatrix_interface/prod
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/solving_strategies/schemes/residualbased_incrementalupdate_static_scheme_slip/ResidualBasedIncrementalUpdateStaticSchemeSlip/Update
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`virtual void RecoverVelocities(ModelPart& rModelPart) const`

# Calls

- [IsSlip](../../../functions/utilities/coordinate_transformation_utilities/IsSlip.md)
- [prod](../../../functions/includes/amatrix_interface/prod.md)

# Called by

- [Update](../../../functions/solving_strategies/schemes/residualbased_incrementalupdate_static_scheme_slip/ResidualBasedIncrementalUpdateStaticSchemeSlip/Update.md)