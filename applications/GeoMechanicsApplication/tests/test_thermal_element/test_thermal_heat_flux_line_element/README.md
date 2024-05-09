# Test Cases for Thermal Point Flux Condition

**Author:** [Mohamed Nabi](https://github.com/mnabideltares)

**Source files:** [Thermal line element with point flux condition](https://github.com/KratosMultiphysics/Kratos/tree/master/applications/GeoMechanicsApplication/tests/test_thermal_element/test_thermal_heat_flux_line_element)

## Case Specification
In this thermal test case, a 1 m high soil column is considered, where the initial temperature in the entire domain is set to 0 $\mathrm{[^\circ C]}$. The top boundary condition is set to be 0 $\mathrm{[^\circ C]}$. A heat flux of 100 $\mathrm{[W/m^2]}$ is set at the bottom boundary. The simulation spans 250 days to allow for a transition from the initial condition to a linear temperature profile between the two sides. This test is conducted for various configurations, including 2D2N, 2D3N, 2D4N, 2D5N, 3D2N and 3D3N line elements. The temperature distribution along the depth is then evaluated with its own result.
The boundary conditions are shown below:

<img src="../documentation_data/test_thermal_point_flux_condition.svg" alt="Visualization of the Boundary conditions" title="Visualization of the Boundary conditions" width="600">

## Results

The picture below illustrates the temperature distribution resulting from the simulation (as an example the 2D3N test is shown below).

<img src="../documentation_data/test_thermal_point_flux_condition_2D3N_result.png" alt="Temperature along the depth at the last time step" title="Temperature along the depth at the last time step" width="600">

These results are associated with the final time step after the solution reaches a steady state. The analytical solution is:

$T = \frac{q}{D} \left( 1 - y \right)$

where $q$ is flux per cross sectional area and $D$ is the dispersion coefficient.

$q = \frac{Q}{A} = \frac{30}{0.3} = 100 \mathrm{[W/m^2]}$ and $D = 1.1328 $ 
$D = n S \lambda^w + (1-n) \lambda^s = 0.36 \times 1 \times 0.48 + (1-0.36) \times 1.5 = 1.1328 \mathrm{[W/m ^\circ C]}$

where, $Q \mathrm{[W]}$ is the total flux, $A \mathrm{[m^2]}$ is the cross section area, $n \mathrm{[-]}$ porosity, $S \mathrm{[-]}$ the saturation, $\lambda^w \mathrm{[W/m ^\circ C]}$ water thermal conductivity and $\lambda^s \mathrm{[W/m ^\circ C]}$ is the soil thermal conductivity.

In this test case, the result at node number 1 at location the bottom of the column, namely at $y = 0 \mathrm{[m]}$, is compared with the analytical solution. The value of the temperature at node 1 is $\frac{100}{1.1328}$ $\mathrm{[^\circ C]}$