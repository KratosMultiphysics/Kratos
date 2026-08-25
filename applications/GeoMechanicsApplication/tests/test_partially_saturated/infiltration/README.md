## Infiltration tests

The geometry and boundary conditions of the infiltration test are shown below. The initial condition is a hydrostatic pressure profile with the reference coordinate at the bottom of the column. This means that at the start, the column contains positive pressures. The boundary condition of $p_w=0$, induces infiltration, which propagates downwards through the column. 

![Geometry and boundary conditions for test_infiltration_pw](test_infiltration_pw/setup_schematic.svg)

Characteristics:
- Zero pressure boundary condition at the top of the model.
- Zero flux boundary condition at the bottom of the model.
- Hydrostatic pressure profile as initial condition with the reference coordinate at y = -2 (at the bottom of the column)
- The elements used are `TransientPwElement2D3N`

### Results
For this problem, linear Pw (water pressure) elements (TransientPwElement2D3N) are used (meaning the displacements are not regarded). The pressure profiles at different time steps can be found in the following image. The red markers depict the asserted pressures, chosen at characteristic positions in the curves.

![Results for Pw elements](test_infiltration_pw/infiltration_from_top_boundary.svg)

### Infiltration tests with 'Staringreeks' materials

Very similar to the set-up described above, the following tests are performed with different materials. The results are shown below.
For each test, the Van Genuchten characteristics are shown in the first image with comparisons to DG-Flow. For the second image, the infiltration from the top boundary is shown with comparisons to DG-Flow, HYDRUS and an external reference solution. The asserted pressures (red markers) are chosen at characteristic positions in the curves.


**Note: one more difference is that the column width of these tests is 0.1m instead of 0.02m**

#### B6
![van_genuchten_characteristics](test_infiltration_pw_caseB6/van_genuchten_characteristics.svg)

![infiltration_from_top_boundary_B6](test_infiltration_pw_caseB6/infiltration_from_top_boundary.svg)


#### B10
![van_genuchten_characteristics](test_infiltration_pw_caseB10/van_genuchten_characteristics.svg)

![infiltration_from_top_boundary_B10](test_infiltration_pw_caseB10/infiltration_from_top_boundary.svg)


#### O4
![van_genuchten_characteristics](test_infiltration_pw_caseO4/van_genuchten_characteristics.svg)

![infiltration_from_top_boundary_O4](test_infiltration_pw_caseO4/infiltration_from_top_boundary.svg)


#### O6
![van_genuchten_characteristics](test_infiltration_pw_caseO6/van_genuchten_characteristics.svg)

![infiltration_from_top_boundary_O6](test_infiltration_pw_caseO6/infiltration_from_top_boundary.svg)

