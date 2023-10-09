## List of modifications in this branch

### Branch information
Branch name: https://github.com/KratosMultiphysics/Kratos/tree/geo/10372-HeatTransfer-temp
Parent branch: master

## Purpose
To add thermal element to GeoMechanicsApplication. This element is decoupled from Pressure and deformation, and hece works independently. Hence, there is no effect of water discharge is added (convection).
Moreover, a condition is added which functions as the normal heat flux condition.

## Steps undertaken
- A new block of thermal elements are added to GeoMechanicsApplications (One element but with different geometries, eg. Triangle2D3N, 6N, etc ...). <br>
applications/GeoMechanicsApplication/geo_mechanics_application.h <br>
Lines 450 to 460

- The thermal elements are registered in GeoMechanicsApplication
applications/GeoMechanicsApplication/geo_mechanics_application.cpp
Lines 211 to 220

- Registered the variables necessary for the thermal element, such as material parameters, time, conditions, etc....
applications/GeoMechanicsApplication/geo_mechanics_application_variables.cpp
Lines 33 to 45

applications/GeoMechanicsApplication/geo_mechanics_application.cpp
lines 324 to 336

applications/GeoMechanicsApplication/geo_mechanics_application_variables.h
lines 60 to 72

- A new element class is added for heat, TransientThermalElement: includes creat, check, intialize, RHS vectors, LHS matrices, etc.
applications/GeoMechanicsApplication/custom_elements/transient_thermal_element.hpp
applications/GeoMechanicsApplication/custom_elements/transient_thermal_element.cpp
Whole files

- The conduction term of the heat equation is multiplied by dispersive matrix. It acts as a constutative law for heat process. So, a constutative law is added, GeoThermalDispersion2DLaw
applications/GeoMechanicsApplication/custom_constitutive/thermal_dispersion_2D_law.hpp
applications/GeoMechanicsApplication/custom_constitutive/thermal_dispersion_2D_law.cpp
Whole files

- The indeces for the dispersive matrix are added (to be consistent with the other laws), INDEX_2D_HEAT_X and INDEX_2D_HEAT_Y
applications/GeoMechanicsApplication/geo_mechanics_application_constants.h
lines 97 to 100

- Normal heat flux condition is added (TNormalFluxCondition) which include RHS, etc... 
applications/GeoMechanicsApplication/custom_conditions/T_normal_flux_condition.hpp
applications/GeoMechanicsApplication/custom_conditions/T_normal_flux_condition.cpp
Whole files

- A thermal condition class is added which include the basic functions of thermal conditions. The thermal conditions are then inherated from this class (TCondition). 
applications/GeoMechanicsApplication/custom_conditions/T_condition.hpp
applications/GeoMechanicsApplication/custom_conditions/T_condition.cpp
Whole files

- A new solver is defined for the thermal element, TSolver, in python. 
applications/GeoMechanicsApplication/python_scripts/geomechanics_T_solver.py
Whole file

- Geomechanics Solver is modified to include the thermal element.
applications/GeoMechanicsApplication/python_scripts/geomechanics_solver.py
Lines 172 to 173, 352 to 357

- GeoMechanics solver warapper is modified to include the thermal element.
applications/GeoMechanicsApplication/python_scripts/geomechanics_solvers_wrapper.py
Lines 24 to 26, Line 30

- A new scheme for the thermal process is added
applications/GeoMechanicsApplication/custom_strategies/schemes/backward_euler_quasistatic_T_scheme.hpp
applications/GeoMechanicsApplication/custom_strategies/schemes/newmark_quasistatic_T_scheme.hpp
Whole files

- Three test cases are added
applications/GeoMechanicsApplication/tests/test_thermal_fixed_temperature_2D3N.gid, to test the thermal application with Dirichlet on the boundary on 3 noded triangles. 
applications/GeoMechanicsApplication/tests/test_thermal_fixed_temperature_2D6N.gid, to test the thermal application with Dirichlet on the boundary on 6 noded triangles. 
applications/GeoMechanicsApplication/tests/test_thermal_heat_flux_2D6N.gid, to test the heat flus boundary on 6 noded triangles (condition uses 3 noded line element). 