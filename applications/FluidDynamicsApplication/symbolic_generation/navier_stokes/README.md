# Navier-Stokes and Stokes elements automatic differentiation

## ELEMENT DESCRIPTION:
Current directory contains the required files for the Automatic Differenctiation (AD) of both the _"symbolic_stokes"_ and _"symbolic_navier_stokes"_ elements of the FluidDynamicsApplication. These elements implement a **Stokes** and a **quasi-incompressible Navier-Stokes** formulation for both **2D** and **3D** cases.

## SYMBOLIC GENERATOR SETTINGS:
*  Dimension to compute: This symbolic generator is valid for both **2D** and **3D** cases. The element has been implemented using a template argument for the problem dimension, so it is advised to set the _dim_to_compute_ flag as "_Both_". In this case the generated .cpp file will contain both **2D** and **3D** implementations.
*  Formulation: As mentioned, the symbolic generator is valid for both Stokes and Navier-Stokes formulations. To switch the formulation set the "_formulation_" variable to either "_Stokes_" or "_NavierStokes_"
*  Linearisation settings: "_FullNR_" considers the convective velocity as _"v-vmesh"_, implying that _v_ is taken into account in the differenctiation of the **LHS** and **RHS**. On the contrary "_Picard_" option (_a.k.a. QuasiNR_) defines the convective velocity as "a", so it is considered as a constant in the differenctiation of the **LHS** and **RHS**. Note that this option is only meaningful in the Navier-Stokes formulation as there is no convective term in the Stokes case.
*  Divide by rho: If _divide_by_rho_ flag is set to "_True_" the mass conservation equation is divided by the density for the sake of having a better conditioned system of equations. It is therefore advised to keep it activated.
*  Artificial compressiblity: If set to "_True_", the time derivative of the density is introduced in the mass conservation equation and rearanged to the pressure time derivative by using the simple state equation $$\partial p/\partial\rho = c^{2}$$, being "_c_" the fluid speed of sound. This adds some extra terms to the usual **Navier-Stokes** equations that are intended to act as a soft artificial compressibility. Such artificial compressibility is controlled by the value of "_c_", meaning that the artificial compressibility terms vanishes as "_c_" tends to infinite. So far this option is only available in the Navier-Stokes formulation. Although the generator is already prepared to include it in the Stokes element, slight modifications would be required in the header and source template files.

## INSTRUCTIONS
Run:
~~~py
python generate_navier_stokes_element.py
~~~
Then, depending on the selected formulation, a file "_symbolic_stokes.cpp_" or "_symbolic_navier_stokes.cpp_" is automatically generated. Such file muest be copied within the "_custom_elements_" folder of the **FluidDynamicsApplication**. The corresponding header file ("symbolic_stokes.h" or "symbolic_navier_stokes.h") is already stored in such folder.
