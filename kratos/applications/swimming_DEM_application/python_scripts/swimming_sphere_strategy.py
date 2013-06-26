import sys
sys.path.insert(0, '../../DEM_application/python_scripts')
sys.path.insert(0, '../../ULFapplication/python_scripts')
from KratosMultiphysics import *
from KratosMultiphysics.IncompressibleFluidApplication import *
from KratosMultiphysics.SwimmingDEMApplication import *
from KratosMultiphysics.PFEMApplication import *

from KratosMultiphysics.StructuralApplication import *
from UD_var import *

# Importing DEM solver strategy

import sphere_strategy as particles

# Importing fluid solver strategy

if (SolverType == "Incompressible_Modified_FracStep"):
    import ulf_frac_swimming as fluid

elif (SolverType == "FracStep"):
    import ulf_frac_swimming as fluid

elif (SolverType == "Quasi_Inc_Constant_Pressure"):
    import ulf_fsi as fluid

elif (SolverType == "Quasi_Inc_Linear_Pressure"):
    import ulf_fsi_inc as fluid

else:
    raise "Solver type not supported: options are fractional_step - \
           modified_frac_steop - quasi_inc_constant_pres - \
           quasi_inc_lin_pres"

def AddVariables(ParticlesModelPart, FluidModelPart):

    particles.AddVariables(ParticlesModelPart)
    fluid.AddVariables(FluidModelPart)

    # HYDRODYNAMICS
    ParticlesModelPart.AddNodalSolutionStepVariable(FLUID_VEL_PROJECTED)
    ParticlesModelPart.AddNodalSolutionStepVariable(FLUID_DENSITY_PROJECTED)
    ParticlesModelPart.AddNodalSolutionStepVariable(PRESSURE_GRAD_PROJECTED)

    # FORCES
    ParticlesModelPart.AddNodalSolutionStepVariable(DRAG_FORCE)
    ParticlesModelPart.AddNodalSolutionStepVariable(BUOYANCY)

def AddDofs(ParticlesModelPart, FluidModelPart):

    particles.AddDofs(ParticlesModelPart)
    fluid.AddDofs(FluidModelPart, compute_reactions)

def DEMStrategy(ParticlesModelPart, domain_size):
    DEM_solver = particles.ExplicitStrategy(ParticlesModelPart, domain_size)

    return DEM_solver

def FluidStrategy(OutFile, FluidOnlyModelPart, FluidModelPart, StructureModelPart, CombinedModelPart, BoxCorner1, BoxCorner2):

    if (SolverType == "Incompressible_Modified_FracStep"):
        fluid_solver = fluid.ULF_FSISolver(OutFile, FluidOnlyModelPart, FluidModelPart, StructureModelPart, CombinedModelPart, FSI, compute_reactions, BoxCorner1, BoxCorner2, domain_size, adaptive_refinement, bulk_modulus, density)

        for node in FluidModelPart.Nodes:
            node.SetSolutionStepValue(BULK_MODULUS, 0, bulk_modulus)
            node.SetSolutionStepValue(DENSITY, 0, density)
            node.SetSolutionStepValue(VISCOSITY, 0, 0.0001)
            node.SetSolutionStepValue(BODY_FORCE_Y, 0, -10.0001)

    if (SolverType == "FracStep"):
        fluid_solver = fluid.ULF_FSISolver(OutFile, FluidOnlyModelPart, FluidModelPart, StructureModelPart, CombinedModelPart, FSI, compute_reactions, BoxCorner1, BoxCorner2, domain_size, adaptive_refinement, bulk_modulus, density)

        for node in FluidModelPart.Nodes:
            node.SetSolutionStepValue(BULK_MODULUS, 0, 0.0)
            node.SetSolutionStepValue(DENSITY, 0, density)

    elif (SolverType == "Quasi_Inc_Constant_Pressure"):
        fluid_solver = fluid.ULF_FSISolver(FluidModelPart, StructureModelPart, CombinedModelPart, compute_reactions, BoxCorner1, BoxCorner2, domain_size, adaptive_refinement)

        for node in FluidModelPart.Nodes:
            node.SetSolutionStepValue(BULK_MODULUS, 0, bulk_modulus)
            node.SetSolutionStepValue(DENSITY, 0, density)

    elif (SolverType == "Quasi_Inc_Linear_Pressure"):
        fluid_solver = fluid.ULF_FSISolver(OutFile, FluidModelPart, StructureModelPart, CombinedModelPart, compute_reactions, BoxCorner1, BoxCorner2, domain_size, adaptive_refinement, bulk_modulus, density)

        for node in FluidModelPart.Nodes:
            node.SetSolutionStepValue(BULK_MODULUS, 0, bulk_modulus)
            node.SetSolutionStepValue(DENSITY, 0, density)

    fluid_solver.alpha_shape = alpha_shape
    fluid_solver.echo_level = 2

    return fluid_solver

class ProjectionModule:

    def __init__(self, FluidModelPart, ParticlesModelPart, DomainSize, NParticlesInDepth):

        self.fluid_model_part     = FluidModelPart
        self.particles_model_part = ParticlesModelPart
        self.n_particles_in_depth = NParticlesInDepth

        if (DomainSize == 3):
            self.projector = BinBasedDEMFluidCoupledMapping3D()
            self.bin_of_objects_fluid = BinBasedFastPointLocator3D(FluidModelPart)

        else:
            self.projector = BinBasedDEMFluidCoupledMapping2D()
            self.bin_of_objects_fluid = BinBasedFastPointLocator2D(FluidModelPart)

    def UpdateDatabase(self, HMin):
        self.bin_of_objects_fluid.UpdateSearchDatabaseAssignedSize(HMin)

    def ProjectFromFluid(self):
        self.projector.InterpolationFromFluidMesh(self.fluid_model_part, self.particles_model_part, self.bin_of_objects_fluid)

    def ProjectFromParticles(self):
        self.projector.InterpolationFromDEMMesh(self.particles_model_part, self.fluid_model_part, self.n_particles_in_depth, self.bin_of_objects_fluid)
