from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics
# Import base class file
from KratosMultiphysics.FluidDynamicsApplication.navier_stokes_ale_fluid_solver import NavierStokesAleFluidSolver
import KratosMultiphysics.ChimeraApplication.python_solvers_wrapper_fluid_chimera as solver_factory

def CreateSolver(model, custom_settings):
    return AleChimeraFluidSolver(model, custom_settings)

class AleChimeraFluidSolver(NavierStokesAleFluidSolver):
    "this class is used only when used defined movement of patch is required"
    def __init__(self,model,project_parameters, parallelism="OpenMP"):
        super(AleChimeraFluidSolver,self).__init__(model,project_parameters, parallelism)

    def _CreateFluidSolver(self, solver_settings, parallelism="OpenMP"):
        return solver_factory.CreateSolverByParameters(self.model, solver_settings, parallelism)

