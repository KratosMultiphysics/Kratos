# Triaxial compression test with 6 noded elements

This test is a drained compression tri-axial test on Mohr-Coulomb model with axi-symmetric 2D6N elements. It mimics a lab test, where soil properties are defined. The test is performed on a cube, consisting of 2 elements, as displayed in the figure below.

![img.png](img.png)

### Setup
The test is performed in two stages, with the following common conditions for both stages:
- Constraints:
  - The displacement in the bottom nodes (5, 7, 9) are fixed in the y direction
  - The displacement in the left nodes (1, 3, 5) are fixed in the x direction.
- Material:
  - The material is described using the Mohr-Coulomb model
- Conditions:
  - An AxisymmetricLineNormalLoadDiffOrderCondition2D3N is added to both the top and right side of the cube (nodes 1, 2, 6 and 6, 8, 9 respectively).

### Stage 1 - Apply a confining stress of -100 kPa, time interval [0, 1]
  - A normal load is applied to the right side of the cube (nodes 6, 8, 9), linearly ramping up from 0 to 100 in the time interval [0, 1].
  - A normal load is applied to the top side of the cube (nodes 1, 2, 6), linearly ramping up from 0 to 100 in the time interval [0, 1].

### Stage 2 - Apply a deviatoric stress of 200 kPa, time interval [1, 1.25]:
  - A constant normal load of 100 is applied to the right side of the cube (nodes 6, 8, 9), during the interval [1.0, 1.25].
  - The displacement of the top nodes (1, 2, 6) is specified in the y direction during the interval [1.0, 2.0] to linearly change from 0.0 to -1.0 meaning at the end-time of 1.25, it will have reached -0.25.

### Checking the results
The calculated effective stresses after the Kratos Geomechanics calculations are compared to the expected solutions:
- After stage 1: The effective stresses in xx, yy and zz are all expected to be -100 in the element integration points.
- After stage 2: The effective stresses in xx and yy are still expected to be -100, while in the zz direction, the expectation is -300 (also here, the comparisons are done in the integration points).





  

    