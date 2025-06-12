# Axisymmetric Navier-Stokes element automatic differentiation

## ELEMENT DESCRIPTION:
Current directory contains the required files for the Automatic Differenctiation (AD) the _"axisymmetric_navier_stokes"_ element of the FluidDynamicsApplication. This element implements an **Newtonian incompressible Navier-Stokes** formulation for the **2D axisymmetric** case.

## AUTOMATIC DIFFERENTIATION SETTINGS:
*  By default, both the linear triangular and linear quadrilateral elements are derived.
*  Divide by rho: If _divide_by_rho_ flag is set to "_True_" the mass conservation equation is divided by the density for the sake of having a better conditioned system of equations. It is therefore advised to keep it activated.
*  ASGS stabilization: The _ASGS_stabilization_ flag is set to "_False_" the subscales are not added to the functional (use only for debugging purposes).
*  Viscous incompressibility error: If _add_incompressibility_error_ is switched to "_True_" the isochoric component of the viscous stress coming from the incompressibility error (i.e. null divergence constraint error) is added to the stress.

## INSTRUCTIONS
Run:
~~~py
python generate_axisymmetric_navier_stokes_element.py
~~~
Then, a file "_axisymmetric_navier_stokes.cpp_" is automatically generated. Such file muest be copied within the "_custom_elements_" folder of the **FluidDynamicsApplication**. The corresponding header file ("axisymmetric_navier_stokes.h") is already stored in such folder.
