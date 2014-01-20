from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# CASE LAUNCHER

from DEM_explicit_solver_var import *

if (ElementType == "SphericParticle3D" or ElementType == "CylinderParticle2D"):
    import spheric_particle_script
elif (ElementType == "SphericContinuumParticle3D" or ElementType == "CylinderContinuumParticle2D"):
    import continuum_spheric_particle_script
