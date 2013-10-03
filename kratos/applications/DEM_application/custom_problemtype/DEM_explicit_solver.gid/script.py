
#CASE LAUNCHER

from DEM_explicit_solver_var import *

if (ElementType == "SphericParticle3D" or ElementType == "CylinderParticle2D"):
    import spheric_particle_script
elif (ElementType == "SphericContinuumParticle3D" or ElementType == "CylinderContinuumParticle2D"):
    import continuum_spheric_particle_script
  
  
