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
  def __init__(self,domain_size,model_part,time_scheme,linear_solver,convergence_criteria,CalculateReactionsFlag,ReformDofSetAtEachStep,MoveMeshFlag,Comm,guess_row_size,periodic_var,dof_per_node=None):
    # add my own builder and solver
    if dof_per_node == None:
        dof_per_node = domain_size+1
    self.builder_and_solver = TrilinosBuilderAndSolverMLPeriodic(Comm,guess_row_size,domain_size,dof_per_node,linear_solver,periodic_var)
    #initialize skiping builder and solver
    SolvingStrategyPython.__init__(self,"NONE",model_part,time_scheme,linear_solver,convergence_criteria,CalculateReactionsFlag,ReformDofSetAtEachStep,MoveMeshFlag,Comm,guess_row_size)
    
