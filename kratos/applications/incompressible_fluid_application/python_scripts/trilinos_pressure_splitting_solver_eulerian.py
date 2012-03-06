try:
 import boost.mpi as mpi
except ImportError:
 import mpi

#importing the Kratos Library
from Kratos import *
from KratosMetisApplication import *
from KratosTrilinosApplication import *
from KratosIncompressibleFluidApplication import *
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

    print "variables for the pressure splitting solver added correctly"
        
def AddDofs(model_part):
    for node in model_part.Nodes:
        #adding dofs
        node.AddDof(VELOCITY_X,REACTION_X);
        node.AddDof(VELOCITY_Y,REACTION_Y);
        node.AddDof(VELOCITY_Z,REACTION_Z);
        node.AddDof(PRESSURE,REACTION_WATER_PRESSURE);
#	node.AddDof(AIR_PRESSURE,REACTION_AIR_PRESSURE);

    mpi.world.barrier()
        
    print "dofs for the pressure splitting solver added correctly"

class PressureSplittingSolver:
    #######################################################################
    def __init__(self,model_part,domain_size):

        self.model_part = model_part
        self.domain_size = domain_size

        self.alpha = -0.3
        self.move_mesh_strategy = 0

        self.Comm = CreateCommunicator()

	#defining the velocity linear solver
        vel_aztec_parameters = ParameterList()
        vel_aztec_parameters.set("AZ_solver","AZ_bicgstab")#"AZ_cg");
        vel_aztec_parameters.set("AZ_kspace",100);
        vel_aztec_parameters.set("AZ_output",32);

##        # preconditioner
##        vel_aztec_parameters.set("AZ_precond","AZ_dom_decomp")
##        vel_aztec_parameters.set("AZ_subdomain_solve","AZ_icc")
##
####        preconditioner_type = "Amesos"
####        preconditioner_parameters = ParameterList()
####        preconditioner_parameters.set("amesos: solver type", "Amesos_Klu");
####        preconditioner_type = "Ifpack_PointRelaxation"
####        preconditioner_parameters = ParameterList()
####        preconditioner_parameters.set("relaxation: type", "Jacobi");
##
####        preconditioner_type = "ILU"
####        preconditioner_parameters = ParameterList()
####        overlap_level = 1
####        nit_max = 1000
####        tol = 1e-6

        # Velocity preconditioner
        vel_preconditioner_type = "ILU"#"IC"
        vel_preconditioner_parameters = ParameterList()
        vel_overlap_level = 1
        vel_nit_max = 1000
        vel_linear_tol = 1e-6

        self.vel_linear_solver =  AztecSolver(\
            vel_aztec_parameters,vel_preconditioner_type,\
            vel_preconditioner_parameters,\
            vel_linear_tol,vel_nit_max,vel_overlap_level);

        # Pressure linear solver
        press_aztec_parameters = ParameterList()
        press_aztec_parameters.set("AZ_solver","AZ_bicgstab")
        press_aztec_parameters.set("AZ_kspace",100);
        press_aztec_parameters.set("AZ_output",32)
        # Pressure preconditioner
        press_preconditioner_type = "ILU"
        press_preconditioner_parameters = ParameterList()
        # Solver parameters
        press_overlap_level = 1
        press_nit_max = 1000
        press_linear_tol = 1e-5

        self.press_linear_solver =  AztecSolver(\
            press_aztec_parameters,press_preconditioner_type,\
            press_preconditioner_parameters,\
            press_linear_tol,press_nit_max,press_overlap_level);
        

        #definition of the convergence criteria
        self.rel_vel_tol = 1e-5
        self.abs_vel_tol = 1e-7
        self.rel_pres_tol = 1e-5
        self.abs_pres_tol = 1e-7

        self.max_iter = 20
                            
        #default settings
        self.echo_level = 1
        self.CalculateReactionFlag = False
        self.ReformDofSetAtEachStep = False
        self.CalculateNormDxFlag = True
        self.MoveMeshFlag = False

        self.dynamic_tau = None
        self.oss_switch  = None
    
        if(domain_size == 2):
            estimate_neighbours = 10
	else:
            estimate_neighbours = 25

        self.guess_row_size = estimate_neighbours * (self.domain_size  + 1)
        self.buildertype="PressureSplitting"
            
##        self.guess_row_size = 15
##        self.buildertype="standard"

        # For Spalart-Allmaras
        self.turbulence_model = None
        self.domain_size = domain_size
            
        
    #######################################################################
    def Initialize(self):

        if self.turbulence_model == None:
            self.time_scheme = TrilinosPredictorCorrectorVelocityBossakScheme\
                               ( self.alpha,self.move_mesh_strategy )
        else:
            self.time_scheme = TrilinosPredictorCorrectorVelocityBossakSchemeTurbulent\
                               (self.alpha,\
                                self.move_mesh_strategy,\
                                self.turbulence_model)

        self.time_scheme.Check(self.model_part)

        # set values for some parameters if provided
        if self.dynamic_tau != None:
            self.model_part.ProcessInfo.SetValue(DYNAMIC_TAU, self.dynamic_tau)
        if self.oss_switch != None:
            self.model_part.ProcessInfo.SetValue(OSS_SWITCH, self.oss_switch )

        #creating the solution strategy
        from trilinos_decoupled_up_strategy_python import DecoupledUPStrategyPython

        self.conv_criteria = TrilinosUPCriteria(self.rel_vel_tol,\
                                                self.abs_vel_tol,\
                                                self.rel_pres_tol,\
                                                self.abs_pres_tol,\
                                                self.Comm)

        
        self.solver = DecoupledUPStrategyPython(\
            self.buildertype,self.model_part,self.time_scheme,\
            self.vel_linear_solver,self.press_linear_solver,self.conv_criteria,\
            self.max_iter,self.CalculateReactionFlag,self.ReformDofSetAtEachStep,\
            self.MoveMeshFlag,self.Comm,self.guess_row_size)

##        self.solver.Check()

##        self.solver = ResidualBasedNewtonRaphsonStrategy(self.model_part,self.time_scheme,self.linear_solver,self.conv_criteria,self.max_iter,self.CalculateReactionFlag, self.ReformDofSetAtEachStep,self.MoveMeshFlag)   
##        (self.solver).SetEchoLevel(self.echo_level)

	                     
    #######################################################################   
    def Solve(self):

        (self.solver).Solve()

       

    #######################################################################   
    def SetEchoLevel(self,level):
        (self.solver).SetEchoLevel(level)
    
    ########################################################################


    def ActivateSmagorinsky(self,c):
        for elem in self.model_part.Elements:
            elem.SetValue(C_SMAGORINSKY,c)

    def ActivateSpalartAllmaras(self,wall_nodes,DES=False):

        number_of_avg_elems = 10
        number_of_avg_nodes = 10
        neighbour_search = FindNodalNeighboursProcess(self.model_part,number_of_avg_elems,number_of_avg_nodes)
        neighbour_search.Execute()

        for node in wall_nodes:
          node.SetValue(IS_VISITED,1.0)
          node.SetSolutionStepValue(DISTANCE,0,0.0)

        # Compute distance function
        if self.domain_size == 2:
            distance_calculator = ParallelDistanceCalculator2D()
        else:
            distance_calculator = ParallelDistanceCalculator3D()

        max_levels = 15
        max_dist = 100.0
        distance_calculator.CalculateDistancesLagrangianSurface\
                            (self.model_part,
                             DISTANCE,
                             NODAL_AREA,
                             max_levels,
                             max_dist)


        import PressureMultiLevelSolver
        turbulence_linear_solver =  PressureMultiLevelSolver.MultilevelLinearSolver(1e-4,1000)

        non_linear_tol = 0.001
        max_it = 3
        reform_dofset = self.ReformDofSetAtEachStep
        time_order = 2

        self.C = CreateCommunicator()
        self.turbulence_model = TrilinosSpalartAllmarasTurbulenceModel(self.C,self.model_part,turbulence_linear_solver,self.domain_size,non_linear_tol,max_it,reform_dofset,time_order)
##        self.turbulence_model = TrilinosSpalartAllmarasTurbulenceModel(self.Comm,self.model_part,turbulence_linear_solver,self.domain_size,non_linear_tol,max_it,reform_dofset,time_order)
        if DES:
            self.turbulence_model.ActivateDES(1.0)




