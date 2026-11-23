---
type: C++ Method
title: ResidualBasedIncrementalAitkenStaticScheme
resource: solving_strategies/schemes/residualbased_incremental_aitken_static_scheme.h#L92-L95
generated:
  by: okf-rs/0.7.0
---

# Signature

`explicit ResidualBasedIncrementalAitkenStaticScheme(Parameters ThisParameters) : ResidualBasedIncrementalAitkenStaticScheme([](Parameters x) -> double {x.ValidateAndAssignDefaults(StaticGetDefaultParameters()); return x["default_omega"].GetDouble(); }(ThisParameters))`