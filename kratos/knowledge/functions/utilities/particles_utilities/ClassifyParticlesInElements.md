---
type: C++ Function
title: ClassifyParticlesInElements
resource: utilities/particles_utilities.h#L144-L191
generated:
  by: okf-rs/0.7.0
---

# Signature

`static void ClassifyParticlesInElements( BinBasedFastPointLocator<TDim>& rLocator, ModelPart& rVolumeModelPart, const ModelPart& rParticlesModelPart, const int NumberOfTypes, const Variable<TScalarType>& rParticleTypeVariable=AUX_INDEX, const Variable<Vector>& rClassificationVectorVariable=MARKER_LABELS, const double SearchTolerance=1e-5 )`