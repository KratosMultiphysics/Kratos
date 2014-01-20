from __future__ import unicode_literals, print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import sys
sys.path.insert(0, '../../DEM_application/python_scripts')
sys.path.insert(0, '../../ULFapplication/python_scripts')

from KratosMultiphysics import *
from KratosMultiphysics.IncompressibleFluidApplication import *
from KratosMultiphysics.SwimmingDEMApplication import *
from KratosMultiphysics.PFEMApplication import *
from KratosMultiphysics.StructuralApplication import *

# Importing DEM solver strategy


class ULFDEMStrategy:

    def __init__(self, creator_destructor, problem_parameters):

        self.param = problem_parameters
        self.creator_destructor = creator_destructor
        self.ImportStrategies()

    def ImportStrategies(self):

        import sphere_strategy as particles

        # Importing fluid solver strategy

        if (self.param.SolverType == "Incompressible_Modified_FracStep"):
            import ulf_frac_swimming as fluid

        elif (self.param.SolverType == "FracStep"):
            import ulf_frac_swimming as fluid

        elif (self.param.SolverType == "Quasi_Inc_Constant_Pressure"):
            import ulf_fsi as fluid

        elif (self.param.SolverType == "Quasi_Inc_Linear_Pressure"):
            import ulf_fsi_inc as fluid

        else:
            raise "Solver type not supported: options are fractional_step - \
                   modified_frac_steop - quasi_inc_constant_pres - \
                   quasi_inc_lin_pres"

        self.DEM = particles
        self.fluid = fluid

    def AddVariables(self, ParticlesModelPart, FluidModelPart):

        self.DEM.AddVariables(ParticlesModelPart, self.param)
        self.fluid.AddVariables(FluidModelPart)

        # HYDRODYNAMICS
        ParticlesModelPart.AddNodalSolutionStepVariable(FLUID_VEL_PROJECTED)
        ParticlesModelPart.AddNodalSolutionStepVariable(FLUID_DENSITY_PROJECTED)
        ParticlesModelPart.AddNodalSolutionStepVariable(PRESSURE_GRAD_PROJECTED)
        ParticlesModelPart.AddNodalSolutionStepVariable(FLUID_VISCOSITY_PROJECTED)

        # FORCES
        ParticlesModelPart.AddNodalSolutionStepVariable(DRAG_FORCE)
        ParticlesModelPart.AddNodalSolutionStepVariable(BUOYANCY)

    def AddDofs(self, ParticlesModelPart, FluidModelPart):

        self.DEM.AddDofs(ParticlesModelPart)
        self.fluid.AddDofs(FluidModelPart, self.param.compute_reactions)

    def DEMStrategy(self, ParticlesModelPart):
        DEM_solver = self.DEM.ExplicitStrategy(ParticlesModelPart, self.particle_destructor, self.param)

        return DEM_solver

    def FluidStrategy(self, OutFile, FluidOnlyModelPart, FluidModelPart, StructureModelPart, CombinedModelPart, BoxCorner1, BoxCorner2):

        if (self.param.SolverType == "Incompressible_Modified_FracStep"):
            fluid_solver = self.fluid.ULF_FSISolver(OutFile, FluidOnlyModelPart, FluidModelPart, StructureModelPart, CombinedModelPart, self.param.FSI, self.param.compute_reactions, BoxCorner1, BoxCorner2, self.param.domain_size, self.param.adaptive_refinement, self.param.bulk_modulus, self.param.density)

            for node in FluidModelPart.Nodes:
                node.SetSolutionStepValue(BULK_MODULUS, 0, self.param.bulk_modulus)
                node.SetSolutionStepValue(DENSITY, 0, self.param.density)
                node.SetSolutionStepValue(VISCOSITY, 0, 0.0001)
                node.SetSolutionStepValue(BODY_FORCE_Y, 0, -10.0001)

        if (self.param.SolverType == "FracStep"):
            fluid_solver = self.fluid.ULF_FSISolver(OutFile, FluidOnlyModelPart, FluidModelPart, StructureModelPart, CombinedModelPart, self.param.FSI, self.param.compute_reactions, BoxCorner1, BoxCorner2, self.param.domain_size, self.param.adaptive_refinement, self.param.bulk_modulus, self.param.density)

            for node in FluidModelPart.Nodes:
                node.SetSolutionStepValue(BULK_MODULUS, 0, 0.0)
                node.SetSolutionStepValue(DENSITY, 0, self.param.density)

        elif (self.param.SolverType == "Quasi_Inc_Constant_Pressure"):
            fluid_solver = self.fluid.ULF_FSISolver(FluidModelPart, StructureModelPart, CombinedModelPart, self.param.compute_reactions, BoxCorner1, BoxCorner2, self.param.domain_size, self.param.adaptive_refinement)

            for node in FluidModelPart.Nodes:
                node.SetSolutionStepValue(BULK_MODULUS, 0, self.param.bulk_modulus)
                node.SetSolutionStepValue(DENSITY, 0, self.param.density)

        elif (self.param.SolverType == "Quasi_Inc_Linear_Pressure"):

            fluid_solver = self.fluid.ULF_FSISolver(OutFile, FluidModelPart, StructureModelPart, CombinedModelPart, self.param.compute_reactions, BoxCorner1, BoxCorner2, self.param.domain_size, self.param.adaptive_refinement, self.param.bulk_modulus, self.param.density)

            for node in FluidModelPart.Nodes:
                node.SetSolutionStepValue(BULK_MODULUS, 0, self.param.bulk_modulus)
                node.SetSolutionStepValue(DENSITY, 0, self.param.density)

        fluid_solver.alpha_shape = self.param.alpha_shape
        fluid_solver.echo_level = 2

        return fluid_solver
