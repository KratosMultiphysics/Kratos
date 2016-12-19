ELEMENT DESCRIPTION:
Current directory contains the documentation for the symbolic derivation of the
"navier_stokes" element. This element includes a formulation of a quasi-incompressible
Navier-Stokes element for both 2D and 3D cases.

SYMBOLIC GENERATOR SETTINGS:
Dimension to compute: This symbolic generator is valid for both 2D and 3D cases.
                      Since the element has been programed with a dimension template in Kratos,
                      it is advised to set the dim_to_compute flag as "Both". In this case the
                      generated .cpp file will contain both 2D and 3D implementations.
Linearisation settings: FullNR considers the convective velocity as "v-vmesh", hence v is taken
                        into account in the derivation of the LHS and RHS. Picard (a.k.a. QuasiNR)
                        considers the convective velocity as "a", thus it is considered as a
                        constant in the derivation of the LHS and RHS.
Artificial compressiblity: If set to true, the time derivative of the density is introduced
                            in the mass conservation equation together with the state equation
                            dp/drho=c^2 (being c the sound velocity). Besides, the velocity
                            divergence is not considered to be 0. These assumptions add some extra
                            terms to the usual Navier-Stokes equations that are intended to act
                            as a soft artificial compressibility, which is controlled by the value
                            of "c". In the derivation of the residual equations the space variations
                            of the density are considered to be close to 0. Besides, the boundary
                            terms arising from the artificial compressibility are also considered to be close to 0.

INSTRUCTIONS
Running "python generate_navier_stokes_element.py" a file "navier_stokes.cpp" is generated
automatically. Such file should be copied within the "custom_elements" folder of the
FluidDynamicsApplication. The corresponding header file "navier_stokes.h", which implements
the element is also stored in the custom_elements folder.
