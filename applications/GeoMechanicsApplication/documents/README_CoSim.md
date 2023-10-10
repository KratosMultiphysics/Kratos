## List of modifications in this branch

### Branch information
Branch name: https://github.com/KratosMultiphysics/Kratos/tree/geo/10372-HeatTransfer-temp-branch-cosim
Parent branch: https://github.com/KratosMultiphysics/Kratos/tree/geo/10372-HeatTransfer-temp

## Purpose
Coupling the thermal element to pressure and deformations. Thuis is done in decoupled-element way. It means that thermal element is kept separetly from the pressure (or pressure-deformation) element, and the coupling is emplemented in the iterative level. Hence a coupling mechanism is needed. Here, the CoSimApp is used.

## Steps undertaken
- A wrapper is added to the CoSimulationApplication. This wrapper cares about the linkes related to the GeoMechechanical python solvers. <br>
Location: applications/CoSimulationApplication/python_scripts/solver_wrappers/kratos/geomechanics_wrapper.py <br>
Whole file

- The discharge is calculated to be used later in convection and dispersion. <br>
applications\GeoMechanicsApplication\custom_elements\transient_thermal_element.cpp <br>
Lines: 666 to 690

- In order to calculate the discharge, from water pressure, the permiability matrix is constructed.
applications\GeoMechanicsApplication\custom_elements\transient_thermal_element.cpp <br>
Lines: 694 to 709

- Convection term is added to the LHS (as matrix) in TransientThermalElement. <br>
applications\GeoMechanicsApplication\custom_elements\transient_thermal_element.cpp <br>
Lines: 372 to 386 <br>
Lines: 505 to 517

- Convection term is added to the RHS (as vector) in TransientThermalElement. <br>
applications\GeoMechanicsApplication\custom_elements\transient_thermal_element.cpp <br>
Lines: 597 to 610 <br>
Lines: 614 to 629

- The required coefficiets are added <br>
applications\GeoMechanicsApplication\custom_elements\transient_thermal_element.cpp <br>
Lines: 654 to 659

- A check for viscosity is added, as viscosity has to be read from thermal material parameters in the case of coupling. <br>
applications\GeoMechanicsApplication\custom_elements\transient_thermal_element.cpp <br>
Lines: 211 to 214

The update for density and viscosity to the thermal element is added <br>
applications\GeoMechanicsApplication\custom_elements\transient_thermal_element.cpp <br>
Lines: 259 to 264

- thermal utilities is added to calculate functions related to thermal effects, here density and viscosity, ThermalUtilities <br>
dapplications\GeoMechanicsApplication\custom_utilities\thermal_utilities.hpp <br>
Whole file

The dispersion matrix is extended to incooperate the effect of water flow <br>
applications\GeoMechanicsApplication\custom_constitutive\thermal_dispersion_2D_law.cpp <br>
Lines: 52 to 63

- The variations in viscosity and density are incoopreted in the pressure element <br>
applications\GeoMechanicsApplication\custom_elements\transient_Pw_element.cpp <br>
Lines: 494 to 497 

- The variation in viscosity and density are incooperated in the following pressure-deformation elements. <br>

  - applications/GeoMechanicsApplication/custom_elements/U_Pw_small_strain_FIC_element.cpp <br>
  Lines: 473 to 477

  - applications/GeoMechanicsApplication/custom_elements/U_Pw_small_strain_element.cpp <br>
  Lines: 1125 to 1129

  - applications/GeoMechanicsApplication/custom_elements/U_Pw_updated_lagrangian_FIC_element.cpp <br>
  Lines: 85 to 89

  - applications/GeoMechanicsApplication/custom_elements/U_Pw_updated_lagrangian_element.cpp <br>
  Lines: 98 to 102

  - applications/GeoMechanicsApplication/custom_elements/small_strain_U_Pw_diff_order_element.cpp <br>
  Lines: 1523 to 1528

  - applications/GeoMechanicsApplication/custom_elements/updated_lagrangian_U_Pw_diff_order_element.cpp <br>
  Lines: 79 to 83

