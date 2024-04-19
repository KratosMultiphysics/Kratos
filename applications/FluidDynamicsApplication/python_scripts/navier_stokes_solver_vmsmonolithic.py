# Importing the Kratos Library
import KratosMultiphysics
from KratosMultiphysics.kratos_utilities import IssueDeprecationWarning

# Import applications
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD

# Import base class file
from KratosMultiphysics.FluidDynamicsApplication.navier_stokes_monolithic_solver import NavierStokesMonolithicSolver

def CreateSolver(model, custom_settings):
    IssueDeprecationWarning('', '"navier_stokes_solver_vmsmonolithic.py" module is deprecated. Please use "navier_stokes_monolithic_solver.py".')
    return NavierStokesMonolithicSolver(model, custom_settings)
