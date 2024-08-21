# Test Cases for Dynamic solution

**Author:** [Mohamed Nabi](https://github.com/mnabideltares)

**Source files:** [Dynamic solution](https://github.com/KratosMultiphysics/Kratos/tree/master/applications/GeoMechanicsApplication/tests/test_dynamic/test_constant_point_load_2d_axi)

## Case Specification
In this test case, a two-dimentional axisymmetric soil block with dimensions of 10 m x 10 m is considered. The pressure is fixed at 0 Pa throughout the entire simulation to eliminate the effect of water pressure.  This is done to allow comparison of our results with published semi-analytical solutions. The bottom side of the domain is fixed, while the sides are allowed to move only in the vertical direction. A sudden force of -1000 N is then applied in the vertical direction at the top surface of the block. The simulation spans 1 second. This test is conducted on a mesh with triangular elements of type 2D3N. The deformation at the surface is then compared with semi-analytical results.

The geometry and boundary conditions are shown below:

<img src="../documentation_data/test_fixed_temperature_boundary_conditions.svg" alt="Visualization of the Boundary conditions" title="Visualization of the Boundary conditions" width="600">

## Results

The picture below illustrates the variation in vertical displacement over time for node 191, which is 3 meters away from the loading point. The results are compared with the semi-analytical solution.

<img src="../documentation_data/test_thermal_fixed_temperature_2D9N_result.png" alt="Temperature within the box width at the last time step" title="Temperature within the box width at the last time step" width="600">

In the semi-analytical solution, there is a singulair point around $\tau = 1.1$. However, this behavior, which is caused by the arrival of the Rayleigh wave, is captured by a peak in the numerical solution. 

Note: To avoid numerical oscillations, the force is gradually applied over a time span of 0.1 seconds. 
