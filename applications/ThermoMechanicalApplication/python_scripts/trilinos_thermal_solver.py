#importing the Kratos Library
from Kratos import *
from KratosThermoMechanicalApplication import *
#from KratosIncompressibleFluidApplication import *
from KratosTrilinosApplication import *

try:
 import boost.mpi as mpi
except ImportError:
 import mpi


def AddVariables(model_part,settings):
    model_part.AddNodalSolutionStepVariable(settings.GetConvectionVariable());
    model_part.AddNodalSolutionStepVariable(settings.GetMeshVelocityVariable());   
    model_part.AddNodalSolutionStepVariable(settings.GetUnknownVariable());    
    model_part.AddNodalSolutionStepVariable(SPECIFIC_HEAT);
    model_part.AddNodalSolutionStepVariable(settings.GetVolumeSourceVariable());
    model_part.AddNodalSolutionStepVariable(settings.GetDensityVariable());
    model_part.AddNodalSolutionStepVariable(settings.GetDiffusionVariable());
    model_part.AddNodalSolutionStepVariable(settings.GetSurfaceSourceVariable());
    model_part.AddNodalSolutionStepVariable(NORMAL);
    
    model_part.AddNodalSolutionStepVariable(PARTITION_INDEX)
    
    if(mpi.rank == 0):
      print "variables for the trilinos thermal solver added correctly"
        
def AddDofs(model_part,settings):
    for node in model_part.Nodes:
	node.AddDof(settings.GetUnknownVariable());
       
    if(mpi.rank == 0):
      print "dofs for the trilinos thermal solver solver added correctly"

class Solver:
    #######################################################################
    def __init__(self,model_part,domain_size,my_settings):

        self.model_part = model_part
        self.domain_size = domain_size
        self.settings = my_settings
        
        self.time_scheme = TrilinosResidualBasedIncrementalUpdateStaticScheme()

        self.Comm = CreateCommunicator()

	self.buildertype="standard"

        #definition of the solvers
        aztec_parameters = ParameterList()
        aztec_parameters.set("AZ_solver","AZ_gmres");
        aztec_parameters.set("AZ_kspace",100);
        aztec_parameters.set("AZ_output","AZ_none");
        preconditioner_type = "ILU"
        preconditioner_parameters = ParameterList()
        overlap_level = 0
        nit_max = 1000
        linear_tol = 1e-9
        self.linear_solver =  AztecSolver(aztec_parameters,preconditioner_type,preconditioner_parameters,linear_tol,nit_max,overlap_level);
        
        #definition of the convergence criteria
        self.conv_criteria = TrilinosDisplacementCriteria(1e-6,1e-9,self.Comm)

        self.CalculateReactionFlag = False
        self.ReformDofSetAtEachStep = False
        self.MoveMeshFlag = False
        self.calculate_norm_dx_flag = False
        self.max_iterations = 1

        if(domain_size == 2):
            self.guess_row_size = 20
        else:
            self.guess_row_size = 45
            

        self.guess_row_size = 18


        
    #######################################################################
    def Initialize(self):
 	(self.model_part.ProcessInfo).SetValue(CONVECTION_DIFFUSION_SETTINGS,self.settings)    
 	
        import trilinos_strategy_python
        self.solver = trilinos_strategy_python.SolvingStrategyPython(self.buildertype,self.model_part,self.time_scheme,self.linear_solver,self.conv_criteria,self.CalculateReactionFlag,self.ReformDofSetAtEachStep,self.MoveMeshFlag,self.Comm,self.guess_row_size)
	                     
    #######################################################################   
    def Solve(self):        
        (self.solver).Solve()      

    #######################################################################   
    def SetEchoLevel(self,level):
        (self.solver).SetEchoLevel(level)
    
    ########################################################################

        




