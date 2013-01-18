#importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.IncompressibleFluidApplication import *
from KratosMultiphysics.FluidDynamicsApplication import *
from KratosMultiphysics.mpi import *
from KratosMultiphysics.MetisApplication import *
from KratosMultiphysics.TrilinosApplication import *
# Check that KratosMultiphysics was imported in the main script
CheckForPreviousImport()


def AddVariables(model_part):
    model_part.AddNodalSolutionStepVariable(VELOCITY);
    model_part.AddNodalSolutionStepVariable(ACCELERATION);
    model_part.AddNodalSolutionStepVariable(MESH_VELOCITY);
    model_part.AddNodalSolutionStepVariable(PRESSURE);
    model_part.AddNodalSolutionStepVariable(IS_STRUCTURE);
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
    model_part.AddNodalSolutionStepVariable(VISCOSITY);
    model_part.AddNodalSolutionStepVariable(DENSITY);
    model_part.AddNodalSolutionStepVariable(BODY_FORCE);
    model_part.AddNodalSolutionStepVariable(NODAL_AREA);
    model_part.AddNodalSolutionStepVariable(NODAL_H);
    model_part.AddNodalSolutionStepVariable(ADVPROJ);
    model_part.AddNodalSolutionStepVariable(DIVPROJ);
    model_part.AddNodalSolutionStepVariable(REACTION); 
    model_part.AddNodalSolutionStepVariable(REACTION_WATER_PRESSURE);
    model_part.AddNodalSolutionStepVariable(EXTERNAL_PRESSURE);
    model_part.AddNodalSolutionStepVariable(FLAG_VARIABLE);
    model_part.AddNodalSolutionStepVariable(NORMAL);
    model_part.AddNodalSolutionStepVariable(Y_WALL);
    model_part.AddNodalSolutionStepVariable(PARTITION_INDEX);
    model_part.AddNodalSolutionStepVariable(MOLECULAR_VISCOSITY)
    model_part.AddNodalSolutionStepVariable(TURBULENT_VISCOSITY)
    model_part.AddNodalSolutionStepVariable(TEMP_CONV_PROJ)
    model_part.AddNodalSolutionStepVariable(DISTANCE)
    model_part.AddNodalSolutionStepVariable(PATCH_INDEX)

    print "variables for the MONOLITHIC_SOLVER_EULERIAN added correctly"
        
def AddDofs(model_part):
    for node in model_part.Nodes:
        #adding dofs
        node.AddDof(VELOCITY_X,REACTION_X);
        node.AddDof(VELOCITY_Y,REACTION_Y);
        node.AddDof(VELOCITY_Z,REACTION_Z);
        node.AddDof(PRESSURE,REACTION_WATER_PRESSURE);
        
    print "dofs for the monolithic solver added correctly"

class MonolithicSolver:
    #######################################################################
    def __init__(self,model_part,domain_size):

        self.model_part = model_part
        self.domain_size = domain_size

        self.alpha = -0.3
        self.move_mesh_strategy = 0

        self.Comm = CreateCommunicator()

        # for Spalart-Allmaras
        self.use_spalart_allmaras = False
        self.wall_nodes = None

        
        self.linear_solver =  TrilinosLinearSolver()
        
        #definition of the convergence criteria
        self.vel_criteria = 1e-3
        self.press_criteria = 1e-3
        self.vel_abs_criteria = 1e-9
        self.press_abs_criteria = 1e-9

        ##self.model_part.ProcessInfo.SetValue(DYNAMIC_TAU, 0.001);

        self.max_iter = 20
                            
        #default settings
        self.echo_level = 1
        self.CalculateReactionFlag = False
        self.ReformDofSetAtEachStep = False
        self.CalculateNormDxFlag = True
        self.MoveMeshFlag = False
    
        if(domain_size == 2):
            estimate_neighbours = 10
            self.guess_row_size = estimate_neighbours * (self.domain_size  + 1)
            #self.buildertype="ML2Dpress"
        else:
            estimate_neighbours = 25
            self.guess_row_size = estimate_neighbours * (self.domain_size  + 1)
            #self.buildertype="ML3Dpress"
            
        #self.guess_row_size = 25
        #self.buildertype="standard"
	#aztec_parameters = ParameterList()
	#aztec_parameters.set("AZ_solver","AZ_gmres");
	#aztec_parameters.set("AZ_kspace",200);
	#aztec_parameters.set("AZ_output","AZ_none");
	#aztec_parameters.set("AZ_output",10);
	#preconditioner_type = "ILU"
	#preconditioner_parameters = ParameterList()
	#preconditioner_parameters.set ("fact: drop tolerance", 1e-9);
	#preconditioner_parameters.set ("fact: level-of-fill", 1);
	#overlap_level = 0
	#nit_max = 1000
	#linear_tol = 1e-9
	#self.linear_solver =  AztecSolver(aztec_parameters,preconditioner_type,preconditioner_parameters,linear_tol,nit_max,overlap_level);

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
      
        ##CHAPUZA to set the non historical value of IS_STRUCTURE correctly... to be improved
        for condition in self.model_part.Conditions:
	    if condition.GetValue(IS_STRUCTURE) == 1.0:
	      for node in condition.GetNodes():
		node.SetSolutionStepValue(IS_STRUCTURE,0,1.0)
	self.model_part.GetCommunicator().AssembleCurrentData(IS_STRUCTURE)
	mpi.world.barrier()
	for node in self.model_part.Nodes:
	  if node.GetSolutionStepValue(IS_STRUCTURE,0) != 0.0:
	    node.SetValue(IS_STRUCTURE,1.0)
	    node.SetSolutionStepValue(IS_STRUCTURE,0,1.0)
        
        ##compute normals "correctly"
        self.normal_calculator = NormalCalculationUtils()
	self.normal_calculator.CalculateOnSimplex(self.model_part,self.domain_size,IS_STRUCTURE)

        # If Spalart-Allmaras: Initialize Spalart-Allmaras solver
        if self.use_spalart_allmaras == True:
            for node in self.wall_nodes:
                node.SetValue(IS_VISITED,1.0)
                node.SetSolutionStepValue(DISTANCE,0,0.0)

            if(self.domain_size == 2):
                self.redistance_utils = ParallelDistanceCalculator2D()
            else:
                self.redistance_utils = ParallelDistanceCalculator3D()

            max_levels = 100
            max_distance = 1000
            self.redistance_utils.CalculateDistances(self.model_part,DISTANCE,NODAL_AREA,max_levels,max_distance)

            non_linear_tol = 0.001
            max_it = 10
            reform_dofset = self.ReformDofSetAtEachStep
            time_order = 2

            turb_aztec_parameters = ParameterList()
            turb_aztec_parameters.set("AZ_solver","AZ_gmres");
            turb_aztec_parameters.set("AZ_kspace",100);
            turb_aztec_parameters.set("AZ_output","AZ_none");

            turb_preconditioner_type = "ILU"
            turb_preconditioner_parameters = ParameterList()
            turb_overlap_level = 0
            turb_nit_max = 1000
            turb_linear_tol = 1e-9

            turb_linear_solver =  AztecSolver(turb_aztec_parameters,turb_preconditioner_type,turb_preconditioner_parameters,turb_linear_tol,turb_nit_max,turb_overlap_level)
            turb_linear_solver.SetScalingType(AztecScalingType.LeftScaling)

            self.turbulence_model = TrilinosSpalartAllmarasTurbulenceModel(self.Comm,self.model_part,turb_linear_solver,self.domain_size,non_linear_tol,max_it,reform_dofset,time_order)
            self.time_scheme = TrilinosPredictorCorrectorVelocityBossakSchemeTurbulent( self.alpha,self.move_mesh_strategy,self.domain_size,self.turbulence_model )
        else: # No turbulence model
            self.time_scheme = TrilinosPredictorCorrectorVelocityBossakSchemeTurbulent( self.alpha,self.move_mesh_strategy,self.domain_size,PATCH_INDEX )

        self.time_scheme.Check(self.model_part)
        
        self.conv_criteria = TrilinosUPCriteria(self.vel_criteria,self.vel_abs_criteria,self.press_criteria,self.press_abs_criteria,self.Comm)

        #creating the solution strategy
        import trilinos_strategy_python_periodic
        self.solver = trilinos_strategy_python_periodic.SolvingStrategyPeriodic(self.domain_size,
                                                                                self.model_part,
                                                                                self.time_scheme,
                                                                                self.linear_solver,
                                                                                self.conv_criteria,
                                                                                self.CalculateReactionFlag,
                                                                                self.ReformDofSetAtEachStep,
                                                                                self.MoveMeshFlag,
                                                                                self.Comm,
                                                                                self.guess_row_size,
                                                                                PATCH_INDEX)
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
    def Clear(self):
        (self.solver).Clear()
        
    ########################################################################
    def ActivateSmagorinsky(self,C):
        for elem in self.model_part.Elements:
            elem.SetValue(C_SMAGORINSKY,C)
        
    ########################################################################
    def ActivateSpalartAllmaras(self,wall_nodes,DES,CDES=1.0):
	self.wall_nodes  = wall_nodes
	self.use_spalart_allmaras = True

