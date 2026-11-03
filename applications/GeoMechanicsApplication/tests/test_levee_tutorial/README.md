# Levee seepage test

## Purpose

This test compares the results of a levee model in Kratos GeoMechanicsApplication with results
from the same model in a different comparable finite-element software package. It
provides a comparison between the two implementations for this particular
idealized setup.

## Model description

The model represents a levee with a clay core. A seepage boundary is prescribed
on the right-hand side of the model, which is not the river-bank side.

The model is deliberately simple:

- It consists of a single analysis stage.
- It has no sheet pile.
- It has no soil-interface elements.
- It is a made-up example and is not a CUR166 case.