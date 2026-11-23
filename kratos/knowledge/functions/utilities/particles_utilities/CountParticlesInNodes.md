---
type: C++ Function
title: CountParticlesInNodes
resource: utilities/particles_utilities.h#L81-L132
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static void CountParticlesInNodes( BinBasedFastPointLocator<TDim>& rLocator, ModelPart& rVolumeModelPart, const ModelPart& rParticlesModelPart, const Variable<double>& rCounterVariable, const double SearchTolerance=1e-5 )`

# Calls

- [GetGeometry](../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)