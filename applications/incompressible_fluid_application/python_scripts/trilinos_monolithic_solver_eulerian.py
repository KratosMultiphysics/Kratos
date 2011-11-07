try:
  import boost.mpi as mpi
except ImportError:
  import mpi

#importing the Kratos Library
from Kratos import *
from KratosMetisApplication import *
from KratosIncompressibleFluidApplication import *
from KratosTrilinosApplication import *
from KratosFluidDynamicsApplication import *
#from KratosStructuralApplication import *


def AddVariables(model_part):
    model_part.AddNodalSolutionStepVariable(VELOCITY);
    model_part.AddNodalSolutionStepVariable(ACCELERATION);
    model_part.AddNodalSolutionStepVariable(MESH_VELOCITY);
    model_part.AddNodalSolutionStepVariable(PRESSURE);
    model_part.AddNodalSolutionStepVariable(AIR_PRESSURE);
    model_part.AddNodalSolutionStepVariable(IS_FLUID);
    model_part.AddNodalSolutionStepVariable(IS_POROUS);
    model_part.AddNodalSolutionStepVariable(IS_STRUCTURE);
    model_part.AddNodalSolutionStepVariable(IS_FREE_SURFACE);
    model_part.AddNodalSolutionStepVariable(IS_INTERFACE);
    model_part.AddNodalSolutionStepVariable(IS_BOUNDARY);
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
    model_part.AddNodalSolutionStepVariable(VISCOSITY);
    model_part.AddNodalSolutionStepVariable(DENSITY);
    model_part.AddNodalSolutionStepVariable(DENSITY_AIR);
    model_part.AddNodalSolutionStepVariable(AIR_SOUND_VELOCITY);
    model_part.AddNodalSolutionStepVariable(SOUND_VELOCITY);
    model_part.AddNodalSolutionStepVariable(BODY_FORCE);
    model_part.AddNodalSolutionStepVariable(NODAL_AREA);
    model_part.AddNodalSolutionStepVariable(NODAL_H);
    model_part.AddNodalSolutionStepVariable(ADVPROJ);
    model_part.AddNodalSolutionStepVariable(DIVPROJ);
    model_part.AddNodalSolutionStepVariable(THAWONE);
    model_part.AddNodalSolutionStepVariable(THAWTWO); 
    model_part.AddNodalSolutionStepVariable(REACTION); 
    model_part.AddNodalSolutionStepVariable(REACTION_WATER_PRESSURE);
    model_part.AddNodalSolutionStepVariable(EXTERNAL_PRESSURE);
    model_part.AddNodalSolutionStepVariable(WATER_PRESSURE);
    model_part.AddNodalSolutionStepVariable(AIR_PRESSURE_DT);
    model_part.AddNodalSolutionStepVariable(ARRHENIUS); 

    model_part.AddNodalSolutionStepVariable(PARTITION_INDEX);
    model_part.AddNodalSolutionStepVariable(NORMAL);

    print "variables for the dynamic structural solution added correctly"
        
def AddDofs(model_part):
    for node in model_part.Nodes:
        #adding dofs
        node.AddDof(VELOCITY_X,REACTION_X);
        node.AddDof(VELOCITY_Y,REACTION_Y);
        node.AddDof(VELOCITY_Z,REACTION_Z);
        node.AddDof(PRESSURE,REACTION_WATER_PRESSURE);
	node.AddDof(AIR_PRESSURE,REACTION_AIR_PRESSURE);

    mpi.world.barrier()
        
    print "dofs for the monolithic solver added correctly"

class MonolithicSolver:
    #######################################################################
    def __init__(self,model_part,domain_size):

        self.model_part = model_part
        self.domain_size = domain_size

        self.alpha = -0.3
        self.move_mesh_strategy = 0

        self.Comm = CreateCommunicator()

        self.turbulence_model = None

        #self.time_scheme = TrilinosPredictorCorrectorVelocityBossakScheme( self.alpha,self.move_mesh_strategy )
        self.linear_solver =  TrilinosLinearSolver()
        
        #definition of the convergence criteria
        self.vel_criteria = 1e-3
        self.press_criteria = 1e-3
        self.vel_abs_criteria = 1e-9
        self.press_abs_criteria = 1e-9

        self.model_part.ProcessInfo.SetValue(DYNAMIC_TAU, 0.001);

        self.max_iter = 20
                            
        #default settings
        self.echo_level = 1
        self.CalculateReactionFlag = False
        self.ReformDofSetAtEachStep = True
        self.CalculateNormDxFlag = True
        self.MoveMeshFlag = False
    
        if(domain_size == 2):
            estimate_neighbours = 10
            self.guess_row_size = estimate_neighbours * (self.domain_size  + 1)
            self.buildertype="ML2Dpress"
        else:
            estimate_neighbours = 25
            self.guess_row_size = estimate_neighbours * (self.domain_size  + 1)
            self.buildertype="ML3Dpress"
            
##        self.guess_row_size = 15
        #self.buildertype="standard"
        #solver_parameters = ParameterList()
        #self.linear_solver =  AmesosSolver("Superludist",solver_parameters);

        ########################################################
        #defining the linear solver
        #self.buildertype="standard"
        #aztec_parameters = ParameterList()
        #aztec_parameters.set("AZ_solver","AZ_gmres");
        #aztec_parameters.set("AZ_kspace",100);
        #aztec_parameters.set("AZ_output",32);

        ##preconditioner_type = "Amesos"
        ##preconditioner_parameters = ParameterList()
        ##preconditioner_parameters.set("amesos: solver type", "Amesos_Klu");

        #preconditioner_type = "ILU"
        #preconditioner_parameters = ParameterList()

        #overlap_level = 0
        #nit_max = 500
        #tol = 1e-6

        #self.linear_solver =  AztecSolver(aztec_parameters,preconditioner_type,preconditioner_parameters,tol,nit_max,overlap_level);
        #self.linear_solver.SetScalingType(AztecScalingType.LeftScaling)
        ##############################################################

        
    #######################################################################
    def Initialize(self):
        self.conv_criteria = TrilinosUPCriteria(self.vel_criteria,self.vel_abs_criteria,self.press_criteria,self.press_abs_criteria,self.Comm)

        # time scheme
        if self.turbulence_model == None:
            self.time_scheme = TrilinosPredictorCorrectorVelocityBossakSchemeTurbulent\
                               ( self.alpha,self.move_mesh_strategy )
        else:
            self.time_scheme = TrilinosPredictorCorrectorVelocityBossakSchemeTurbulent\
                               (self.alpha,\
                                self.move_mesh_strategy,\
                                self.turbulence_model)

        #creating the solution strategy
        import trilinos_strategy_python
        self.solver = trilinos_strategy_python.SolvingStrategyPython(self.buildertype,self.model_part,self.time_scheme,self.linear_solver,self.conv_criteria,self.CalculateReactionFlag,self.ReformDofSetAtEachStep,self.MoveMeshFlag,self.Comm,self.guess_row_size)
        self.solver.max_iter = self.max_iter

##        self.solver = ResidualBasedNewtonRaphsonStrategy(self.model_part,self.time_scheme,self.linear_solver,self.conv_criteria,self.max_iter,self.CalculateReactionFlag, self.ReformDofSetAtEachStep,self.MoveMeshFlag)   
        (self.solver).SetEchoLevel(self.echo_level)

	                     
    #######################################################################   
    def Solve(self):

        (self.solver).Solve()

       

    #######################################################################   
    def SetEchoLevel(self,level):
        (self.solver).SetEchoLevel(level)
    
    ########################################################################

        




