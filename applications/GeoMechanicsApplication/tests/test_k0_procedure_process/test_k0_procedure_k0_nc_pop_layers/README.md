# Test K<sub>0</sub> procedure with horizontal layers, K_0 normal consolidation with POP using layers

**Author:** [Jonathan Nuttall](https://github.com/mcgicjn2)

## Case Specification
A test for the K<sub>0</sub> procedure process.  The plane strain model consists of three horizontal soil layers.  The phreatic line coincides with the top boundary of the top soil layer.  Each soil layer has a linear elastic constitutive law associated with it.

K<sub>0</sub> is calculated using K<sub>0</sub><sup>nc</sup> and the POP method.  This also corrects for $\nu_{ur}.

The aim of the test is to verify that the computed effective stress distribution matches the one obtained with comparative software.
