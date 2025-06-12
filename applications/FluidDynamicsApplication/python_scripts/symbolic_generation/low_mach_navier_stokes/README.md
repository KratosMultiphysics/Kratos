# Low Mach Navier-Stokes element automatic differentiation

## ELEMENT DESCRIPTION:
Current directory contains the required files for the Automatic Differenctiation (AD) the _"low_mach_navier_stokes"_ element of the FluidDynamicsApplication. This element implements an **low Mach approximation Navier-Stokes** formulation for the **2D** and **3D** cases.

## AUTOMATIC DIFFERENTIATION SETTINGS:
*  By default, both the linear triangular and linear quadrilateral elements as well as linear tetrahedron and hexahedra are derived.
*  ASGS stabilization: The _ASGS_stabilization_ flag is set to "_False_" the subscales are not added to the functional (use only for debugging purposes).

## INSTRUCTIONS
Run:
~~~py
python generate_low_mach_navier_stokes_element.py
~~~
Then, a file "_low_mach_navier_stokes.cpp_" is automatically generated. Such file muest be copied within the "_custom_elements_" folder of the **FluidDynamicsApplication**. The corresponding header file ("low_mach_navier_stokes.h") is already stored in such folder.
