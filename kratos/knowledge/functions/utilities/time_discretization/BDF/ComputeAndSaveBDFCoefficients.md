---
type: C++ Method
title: ComputeAndSaveBDFCoefficients
resource: utilities/time_discretization.cpp#L58-L78
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/solving_strategies/schemes/residual_based_bdf_scheme/ResidualBasedBDFScheme/InitializeSolutionStep
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void BDF::ComputeAndSaveBDFCoefficients(ProcessInfo &rProcessInfo) const`

# Called by

- [InitializeSolutionStep](../../../../functions/solving_strategies/schemes/residual_based_bdf_scheme/ResidualBasedBDFScheme/InitializeSolutionStep.md)