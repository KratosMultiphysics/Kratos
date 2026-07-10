Current directory contains the documentation for the symbolic derivation of the 
"stokes_3D_twofluid" element, which includes a formulation of a stokes element
including a pressure enrichment to allow dealing with stokes + air

running 
python generate_stokes_twofluid_element.py

a file 
"stokes_3D_twofluid.cpp" is generated automatically
such file should be copied within the "custom_elements" folder

the corresponding header file "stokes_3D_twofluid.h", which implements the element subdivision and the local static condensation, is contained in the 
is also stored in the custom_elements folder

