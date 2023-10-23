## List of modifications in this branch

### Branch information
Branch name: https://github.com/KratosMultiphysics/Kratos/tree/geo/10372-HeatTransfer-temp-branch-cosim
Parent branch: https://github.com/KratosMultiphysics/Kratos/tree/geo/10372-HeatTransfer-temp

## Purpose
Adding an advanced boundary condition, namely micro-climate to the thermal part. This is done bt means of adding new condition where the right and lefthand sides of the conditions are calculated.

## Steps undertaken
- Micro-climate condition is registered in the GeomechanicsApplication. <br>
Location: applications\GeoMechanicsApplication\geo_mechanics_application.h
Lines 499 to 507
applications\GeoMechanicsApplication\geo_mechanics_application.cpp
Lines 254 to 262

- List of variables necessary for micro-climate boundary is registerd in the GeoMechanicsApplication. <br>
applications/GeoMechanicsApplication\geo_mechanics_application.cpp <br>
Lines 356 to 367 <br>
applications\GeoMechanicsApplication\geo_mechanics_application_variables.cpp <br>
46 to 57 <br>
applications\GeoMechanicsApplication\custom_python\geo_mechanics_python_application.cpp <br>
61 to 65 <br>
d:\Kratos\10372-HeatTransfer-temp-branch-climate\applications\GeoMechanicsApplication\python_scripts\geomechanics_solver.py <br>
358 to 362

- A new condition is defined, as class. <br>
applications/GeoMechanicsApplication/custom_conditions/T_microclimate_flux_condition.hpp <br>
applications/GeoMechanicsApplication/custom_conditions/T_normal_flux_condition.cpp <br>
Whole files