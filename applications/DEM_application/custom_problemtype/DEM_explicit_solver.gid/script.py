from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# CASE LAUNCHER

from DEM_explicit_solver_var import *

print (ElementType)
if (ElementType == "SphericPartDEMElement3D" or ElementType == "CylinderPartDEMElement2D"):
    import spheric_particle_script
    
    
elif (ElementType == "SphericContPartDEMElement3D" or ElementType == "CylinderContPartDEMElement2D"):
    import continuum_spheric_particle_script
    
