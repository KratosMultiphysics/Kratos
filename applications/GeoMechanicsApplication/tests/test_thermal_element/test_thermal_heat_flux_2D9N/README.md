# Test Cases for Thermal heat flux

**Author:** [Marjan Fathian](https://github.com/indigocoral)

**Source files:** [Thermal element tests](https://github.com/KratosMultiphysics/Kratos/tree/master/applications/GeoMechanicsApplication/tests/test_thermal_element)


## Case Specification
In this thermal heat flux test case, a 1 x 1 m soil box is considered, with a heat flux applied as a condition to the bottom side and the top temperature set to 0 degrees. The simulation spans 100 days to allow for a transition from an exponential to a linear temperature profile between the two sides. This test is conducted for various configurations, including 2D3N, 2D4N, 2D6N, 2D8N, 2D9N, 2D10N, 2D15N, 3D4N, and 3D10N. The temperature distribution within the box is then assessed using an analytical solution.
The boundary conditions are shown below:

<img src="data/test_heat_flux_boundary_conditions.png" width="300">

_Visualization of the Boundary conditions[1]_


## Results

The animated figure below illustrates the temperature contours resulting from the simulation (as an example the 2D9N test is shown below).

<img src="data/test_thermal_heat_flux_2D9N_result.gif" width="500">

_Temperature contours[2]_


