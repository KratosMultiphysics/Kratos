# Importing the Kratos Library
from KratosMultiphysics.kratos_utilities import IssueDeprecationWarning

# Import current class file
from KratosMultiphysics.FluidDynamicsApplication.trilinos_navier_stokes_monolithic_solver import TrilinosNavierStokesMonolithicSolver

def CreateSolver(model, custom_settings):
    IssueDeprecationWarning('', '"trilinos_navier_stokes_solver_vmsmonolithic.py" module is deprecated. Please use "trilinos_navier_stokes_monolithic_solver.py".')
    return TrilinosNavierStokesMonolithicSolver(model, custom_settings)