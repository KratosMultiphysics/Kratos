# Test Cases for Thermal Line Element temperature

**Author:** [Mohamed Nabi](https://github.com/mnabideltares)

**Source files:** [Thermal line element with fixed temperature](https://github.com/KratosMultiphysics/Kratos/tree/master/applications/GeoMechanicsApplication/tests/test_thermal_element/test_thermal_line_element)

## Case Specification
In this thermal test case, a 3 m deep soil is considered, with everywhere set to 0 $\mathrm{[^\circ C]}$ then a sudden jump at the top and bottom boundary temperatures are given. These are -10 $\mathrm{[^\circ C]}$ at bottom and 100 $\mathrm{[^\circ C]}$ at the top boundaries. The simulation spans 1000 days to allow for a transition from an exponential to a linear temperature profile between the two sides. This test is conducted for various configurations, including 2D2N, 2D3N, 2D4N, 2D5N, 3D2N and 3D3N line elements. The temperature distribution along the depth is then evaluated with its own result.
The boundary conditions are shown below:

<img src="../documentation_data/test_thermal_line_element.svg" alt="Visualization of the Boundary conditions" title="Visualization of the Boundary conditions" width="600">

## Results

The picture below illustrates the temperature contours resulting from the simulation (as an example the 2D3N test is shown below).

<img src="../documentation_data/test_thermal_line_element_2D3N_result.png" alt="Temperature along the depth at the last time step" title="Temperature along the depth at the last time step" width="600">

These results are associated with the final time step after the solution reaches a steady state. The analytical solution is:

$T = \frac{110}{3} y + 100$

In this test case, the result at node number 3 at location $y = -2$ is compared with the analytical solution. The value of the temperature at node 3 is 26.6666666... $\mathrm{[^\circ C]}$