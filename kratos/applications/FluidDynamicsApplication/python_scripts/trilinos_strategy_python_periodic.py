#importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.mpi import *
from KratosMultiphysics.TrilinosApplication import *
from KratosMultiphysics.FluidDynamicsApplication import *
# Check that KratosMultiphysics was imported in the main script
CheckForPreviousImport()

# import base class
from trilinos_strategy_python import SolvingStrategyPython

class SolvingStrategyPeriodic(SolvingStrategyPython):
  def __init__(self,domain_size,model_part,time_scheme,linear_solver,convergence_criteria,CalculateReactionsFlag,ReformDofSetAtEachStep,MoveMeshFlag,Comm,guess_row_size,periodic_var):
    # add my own builder and solver
    self.builder_and_solver = TrilinosBuilderAndSolverMLPeriodic(Comm,guess_row_size,domain_size,linear_solver,periodic_var)
    #initialize skiping builder and solver
    SolvingStrategyPython.__init__(self,"NONE",model_part,time_scheme,linear_solver,convergence_criteria,CalculateReactionsFlag,ReformDofSetAtEachStep,MoveMeshFlag,Comm,guess_row_size)
    
