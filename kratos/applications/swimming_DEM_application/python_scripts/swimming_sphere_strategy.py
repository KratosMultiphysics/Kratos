import sys
from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *
sys.path.insert(0, '../../DEM_application/python_scripts')
sys.path.insert(0, '../../ULFapplication/python_scripts')
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
    fluid_solver.Initialize()
    fluid_solver.echo_level = 2

    return fluid_solver
#def FluidStrategy(FluidModelPart, DomainSize, SolverType):
