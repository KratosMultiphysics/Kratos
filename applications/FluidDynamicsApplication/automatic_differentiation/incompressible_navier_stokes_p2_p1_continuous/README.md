# Incompressible Navier-Stokes P2/P1 continuous automatic differentiation

## ELEMENT DESCRIPTION:
Current directory contains the required files for the Automatic Differenctiation (AD) of the _"incompressible_navier_stokes_p2_p1_continuous"_ element of the FluidDynamicsApplication. This element implements an **incompressible Navier-Stokes** formulation with a P2/P1 continuous (Taylor-Hood) inf-sup stable interpolation pair for both **2D** and **3D** cases. Quasi-static Algebraic Sub-Grid Scales (ASGS) are used for the stabilization of the convection.

## SYMBOLIC GENERATOR SETTINGS:
*  Dimension to compute: this symbolic generator is valid for both **2D** and **3D** cases. The element has been implemented using a template argument for the problem dimension, so it is advised to set the _dim_to_compute_ flag as "_Both_". In this case the generated .cpp file will contain both **2D** and **3D** implementations.
*  Linearisation settings: "_FullNR_" considers the convective velocity as _"v-vmesh"_, implying that _v_ is taken into account in the differenctiation of the **LHS** and **RHS**. On the contrary "_Picard_" option (_a.k.a. QuasiNR_) defines the convective velocity as "a", so it is considered as a constant in the differenctiation of the **LHS** and **RHS**.
*  Add pressure subscale: if set to _True_ the pressure subscale component is added to the momentum equation to effectively get the div-div stabilization term. Though it is not required, this term is known to greatly improve the solution in presence of inf-sup stable approximants (see a detailed discussion in [here](https://doi.org/10.1016/j.cma.2016.02.026)).

## INSTRUCTIONS
Run:
~~~py
python generate_incompressible_navier_stokes_p2_p1_continuous_element.py
~~~
Then, a file "_incompressible_navier_stokes_p2_p1_continuous.cpp_" is automatically generated. Such file muest be moved to the "_custom_elements_" folder of the **FluidDynamicsApplication**. The corresponding header file ("incompressible_navier_stokes_p2_p1_continuous.h") is already stored in there.
