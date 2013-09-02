# -*- coding: utf-8 -*-
from KratosMultiphysics import *
from KratosMultiphysics.mpi import *
from KratosMultiphysics.IncompressibleFluidApplication import *
from KratosMultiphysics.FluidDynamicsApplication import *
from KratosMultiphysics.MetisApplication import *
from KratosMultiphysics.TrilinosApplication import *
# Check that KratosMultiphysics was imported in the main script
CheckForPreviousImport()

def AddVariables(model_part, config=None):
    model_part.AddNodalSolutionStepVariable(VELOCITY);
    model_part.AddNodalSolutionStepVariable(FRACT_VEL);
    model_part.AddNodalSolutionStepVariable(MESH_VELOCITY);
    model_part.AddNodalSolutionStepVariable(PRESSURE);
    model_part.AddNodalSolutionStepVariable(PRESSURE_OLD_IT);
    model_part.AddNodalSolutionStepVariable(PRESS_PROJ);
    model_part.AddNodalSolutionStepVariable(CONV_PROJ);
    model_part.AddNodalSolutionStepVariable(DIVPROJ);
    model_part.AddNodalSolutionStepVariable(NODAL_AREA);
    model_part.AddNodalSolutionStepVariable(BODY_FORCE);
    model_part.AddNodalSolutionStepVariable(DENSITY);
    model_part.AddNodalSolutionStepVariable(VISCOSITY);
    model_part.AddNodalSolutionStepVariable(FLAG_VARIABLE);
    model_part.AddNodalSolutionStepVariable(IS_STRUCTURE);
    model_part.AddNodalSolutionStepVariable(PARTITION_INDEX)
    model_part.AddNodalSolutionStepVariable(REACTION)
    model_part.AddNodalSolutionStepVariable(NORMAL)
    model_part.AddNodalSolutionStepVariable(Y_WALL);
    model_part.AddNodalSolutionStepVariable(MOLECULAR_VISCOSITY)
    model_part.AddNodalSolutionStepVariable(TURBULENT_VISCOSITY)
    model_part.AddNodalSolutionStepVariable(TEMP_CONV_PROJ)
    model_part.AddNodalSolutionStepVariable(DISTANCE)
    mpi.world.barrier()
    if mpi.rank == 0:
        print "variables for the trilinos fractional step solver added correctly"

def AddDofs(model_part, config=None):
  
    for node in model_part.Nodes:
        #adding dofs
        node.AddDof(PRESSURE);
        node.AddDof(VELOCITY_X,REACTION_X);
        node.AddDof(VELOCITY_Y,REACTION_Y);
        node.AddDof(VELOCITY_Z,REACTION_Z);
    mpi.world.barrier()

    if mpi.rank == 0:
        print "dofs for the trilinos fractional step solver added correctly"


class IncompressibleFluidSolver:
    
    def __init__(self,model_part,domain_size):

        #neighbour search
        number_of_avg_elems = 10
        number_of_avg_nodes = 10
        self.neighbour_search = FindNodalNeighboursProcess(model_part,number_of_avg_elems,number_of_avg_nodes)

        self.model_part = model_part
        self.domain_size = domain_size

        #assignation of parameters to be used
        self.vel_toll = float(0.001);
        self.press_toll = float(0.001);
        self.max_vel_its = 6;
        self.max_press_its = 3;
        self.time_order = 2;
        self.compute_reactions = False;
        self.ReformDofAtEachIteration = False; 
        self.CalculateNormDxFlag = True;
        self.predictor_corrector = False;
        self.use_dt_in_stabilization = False

        self.echo_level = 0

        self.step = 1
        self.projections_are_initialized = False;

        self.dynamic_tau = 0.001
        self.activate_tau2 = False
        
        self.Comm = CreateCommunicator()


##        ########################################################
        #defining the linear solver
        velocity_aztec_parameters = ParameterList()
        velocity_aztec_parameters.set("AZ_solver","AZ_bicgstab");
        velocity_aztec_parameters.set("AZ_kspace",100);
        velocity_aztec_parameters.set("AZ_output","AZ_none");
       #velocity_aztec_parameters.set("AZ_output",32);

##        preconditioner_type = "Amesos"
##        preconditioner_parameters = ParameterList()
##        preconditioner_parameters.set("amesos: solver type", "Amesos_Klu");
##        preconditioner_type = "Ifpack_PointRelaxation"
##        preconditioner_parameters = ParameterList()
##        preconditioner_parameters.set("relaxation: type", "Jacobi");

##        preconditioner_type = "ILU"
##        preconditioner_parameters = ParameterList()
##        overlap_level = 1
##        nit_max = 1000
##        tol = 1e-6

        velocity_preconditioner_type = "ILU"
        velocity_preconditioner_parameters = ParameterList()
        velocity_overlap_level = 0
        velocity_nit_max = 1000
        velocity_linear_tol = 1e-6
       

        
        self.velocity_linear_solver =  AztecSolver(velocity_aztec_parameters,velocity_preconditioner_type,velocity_preconditioner_parameters,velocity_linear_tol,velocity_nit_max,velocity_overlap_level);
        #self.velocity_linear_solver.SetScalingType(AztecScalingType.NoScaling)
        #self.pressure_linear_solver.SetScalingType(AztecScalingType.NoScaling) 
        self.velocity_linear_solver.SetScalingType(AztecScalingType.LeftScaling)
        #self.velocity_linear_solver.SetScalingType(AztecScalingType.SymmetricScaling)
        #self.pressure_linear_solver.SetScalingType(AztecScalingType.SymmetricScaling) 
        

        pressure_nit_max = 1000
        pressure_linear_tol = 1e-6
        
        #pressure_aztec_parameters = ParameterList()
        #pressure_aztec_parameters.set("AZ_solver","AZ_bicgstab");
        ##pressure_preconditioner_type = "IC"
        #pressure_preconditioner_type = "ILU"
        ##pressure_preconditioner_type = "AZ_none"
        #pressure_aztec_parameters.set("AZ_output",32);
        ##pressure_aztec_parameters.set("AZ_output","AZ_none");
        #pressure_preconditioner_parameters = ParameterList()
        #pressure_overlap_level = 0
        #self.pressure_linear_solver =  AztecSolver(pressure_aztec_parameters,pressure_preconditioner_type,pressure_preconditioner_parameters,pressure_linear_tol,pressure_nit_max,pressure_overlap_level);
        
        import PressureMultiLevelSolver
        self.pressure_linear_solver =  PressureMultiLevelSolver.MultilevelLinearSolver(pressure_linear_tol,pressure_nit_max)
        #self.pressure_linear_solver.SetScalingType(AztecScalingType.LeftScaling) 
        
        
        
        ##handling slip condition
##        self.slip_conditions_initialized = False
##        self.create_slip_conditions = GenerateSlipConditionProcess(self.model_part,domain_size)

        self.use_slip_conditions = False
        self.use_spalart_allmaras = False
        self.use_des = False

##        ##############################################################

##        vel_solver_parameters = ParameterList()
##        self.velocity_linear_solver =  AmesosSolver("Superludist",vel_solver_parameters);
##
##        press_solver_parameters = ParameterList()
##        self.pressure_linear_solver =  AmesosSolver("Superludist",press_solver_parameters);



    def Initialize(self):
        (self.neighbour_search).Execute()
        
        self.model_part.ProcessInfo.SetValue(DYNAMIC_TAU, self.dynamic_tau);
        self.model_part.ProcessInfo.SetValue(ACTIVATE_TAU2, self.activate_tau2);

        # check if slip conditions are defined
        slip_cond_count = 0
        if self.use_slip_conditions == False:
            for cond in self.model_part.Conditions:
                if cond.GetValue(IS_STRUCTURE) != 0.0:
                    slip_cond_count += 1
                    break
        all_results = mpi.allgather(mpi.world,slip_cond_count)
        for count in all_results:
            slip_cond_count += count
        if slip_cond_count > 0:
            self.use_slip_conditions = True

        # if we use slip conditions, calculate normals on the boundary
        if self.use_slip_conditions == True:
            self.normal_util = NormalCalculationUtils()
            self.normal_util.CalculateOnSimplex(self.model_part,self.domain_size,IS_STRUCTURE)
        
#        self.solver = ResidualBasedFluidStrategyCoupled(self.model_part,self.velocity_linear_solver,self.pressure_linear_solver,self.compute_reactions,self.ReformDofAtEachIteration,self.CalculateNormDxFlag,self.vel_toll,self.press_toll,self.max_vel_its,self.max_press_its, self.time_order,self.domain_size, self.laplacian_form, self.predictor_corrector)   
#        print "in python: okkio using Coupled Strategy"
#        self.solver = ResidualBasedFluidStrategy(self.model_part,self.velocity_linear_solver,self.pressure_linear_solver,self.compute_reactions,self.ReformDofAtEachIteration,self.CalculateNormDxFlag,self.vel_toll,self.press_toll,self.max_vel_its,self.max_press_its, self.time_order,self.domain_size, self.laplacian_form, self.predictor_corrector)

        
##        print self.Comm

        MoveMeshFlag = False
##        solver_configuration = TrilinosFractionalStepConfiguration(self.Comm,self.model_part,self.velocity_linear_solver,self.pressure_linear_solver,self.domain_size,self.laplacian_form,self.use_dt_in_stabilization )

##        print "solver configuration created correctly"
##        self.solver = TrilinosFractionalStepStrategy( self.model_part, solver_configuration, self.ReformDofAtEachIteration, self.vel_toll, self.press_toll, self.max_vel_its, self.max_press_its, self.time_order, self.domain_size,self.predictor_corrector)

##        self.solver = TrilinosFSStrategy(self.model_part,
##                                         self.velocity_linear_solver,
##                                         self.pressure_linear_solver,
##                                         MoveMeshFlag,
##                                         self.ReformDofAtEachIteration,
##                                         self.vel_toll,
##                                         self.press_toll,
##                                         self.max_vel_its,
##                                         self.max_press_its,
##                                         self.time_order,
##                                         self.domain_size,
##                                         self.predictor_corrector)

        self.solver_settings = TrilinosFractionalStepSettings(self.Comm,
                                              self.model_part,
                                              self.domain_size,
                                              self.time_order,
                                              self.use_slip_conditions,
                                              MoveMeshFlag,
                                              self.ReformDofAtEachIteration)

        self.solver_settings.SetStrategy(TrilinosStrategyLabel.Velocity,
                                         self.velocity_linear_solver,
                                         self.vel_toll,
                                         self.max_vel_its)

        self.solver_settings.SetStrategy(TrilinosStrategyLabel.Pressure,
                                         self.pressure_linear_solver,
                                         self.press_toll,
                                         self.max_press_its)

        if self.use_spalart_allmaras:
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
            reform_dofset = self.ReformDofAtEachIteration
            time_order = self.time_order

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

            self.solver_settings.SetTurbulenceModel(
                TrilinosTurbulenceModelLabel.SpalartAllmaras,
                turb_linear_solver,
                non_linear_tol,
                max_it)

            self.solver_settings.SetEchoLevel(self.echo_level)

        self.solver = TrilinosFSStrategy(self.model_part,
                                         self.solver_settings,
                                         self.predictor_corrector)


##        ##generating the slip conditions
##        self.create_slip_conditions.Execute()
##        #(self.solver).SetSlipProcess(self.create_slip_conditions);
##        self.slip_conditions_initialized = True

##        (self.solver).SetEchoLevel(self.echo_level)



##        (self.solver).SetEchoLevel(self.echo_level)
        print "finished initialization of the fluid strategy"
        
   
    def Solve(self):
        if(self.ReformDofAtEachIteration == True):
            self.slip_conditions_initialized = False
            (self.neighbour_search).Execute()
##        if(self.slip_conditions_initialized == False):
##            self.create_slip_conditions.Execute()
##            #(self.solver).SetSlipProcess(self.create_slip_conditions);
##            self.slip_conditions_initialized = True

        (self.solver).Solve()

        if(self.compute_reactions == True):
            self.solver.CalculateReactions()
        #(self.solver).ApplyFractionalVelocityFixity()
        #(self.solver).InitializeFractionalStep(self.step, self.time_order);
        #(self.solver).InitializeProjections(self.step,self.projections_are_initialized);
        #self.projections_are_initialized = True;

        #(self.solver).AssignInitialStepValues();

        #(self.solver).SolveStep1(self.vel_toll,self.max_vel_its)
        #(self.solver).SolveStep2();
        #(self.solver).ActOnLonelyNodes();
        #(self.solver).SolveStep3();
        #(self.solver).SolveStep4();

    def Clear(self):
        (self.solver).Clear()
        self.slip_conditions_initialized = True

    def WriteRestartFile(self,FileName):
        backupfile = open(FileName+".py",'w')
        backupfile.write( "from KratosMultiphysics import *\n");
        backupfile.write( "def Restart(NODES):\n" )
        
        import restart_utilities
        restart_utilities.PrintRestart_ScalarVariable_PyFormat(VELOCITY_X,"VELOCITY_X",self.model_part.Nodes,backupfile)
        restart_utilities.PrintRestart_ScalarVariable_PyFormat(VELOCITY_Y,"VELOCITY_Y",self.model_part.Nodes,backupfile)
        restart_utilities.PrintRestart_ScalarVariable_PyFormat(VELOCITY_Z,"VELOCITY_Z",self.model_part.Nodes,backupfile)
        restart_utilities.PrintRestart_ScalarVariable_PyFormat(PRESSURE,"PRESSURE",self.model_part.Nodes,backupfile)
        restart_utilities.PrintRestart_ScalarVariable_PyFormat(DENSITY,"DENSITY",self.model_part.Nodes,backupfile)
        restart_utilities.PrintRestart_ScalarVariable_PyFormat(VISCOSITY,"VISCOSITY",self.model_part.Nodes,backupfile)

        restart_utilities.PrintRestartFixity_PyFormat(VELOCITY_X,"VELOCITY_X",self.model_part.Nodes,backupfile)
        restart_utilities.PrintRestartFixity_PyFormat(VELOCITY_Y,"VELOCITY_Y",self.model_part.Nodes,backupfile)
        restart_utilities.PrintRestartFixity_PyFormat(VELOCITY_Z,"VELOCITY_Z",self.model_part.Nodes,backupfile)
        restart_utilities.PrintRestartFixity_PyFormat(PRESSURE,"PRESSURE",self.model_part.Nodes,backupfile)
        
        backupfile.close()


    def ActivateSpalartAllmaras(self,wall_nodes,DES,CDES=1.0):
        self.wall_nodes  = wall_nodes
        self.use_spalart_allmaras = True
        #for node in wall_nodes:
            #node.SetValue(IS_VISITED,1.0)

        #if(self.domain_size == 2):
            #self.redistance_utils = ParallelDistanceCalculator2D()
        #else:
            #self.redistance_utils = ParallelDistanceCalculator3D()

        #max_levels = 100
        #max_distance = 1000
        #self.redistance_utils.CalculateDistances(self.model_part,DISTANCE,NODAL_AREA,max_levels,max_distance)

        #non_linear_tol = 0.001
        #max_it = 10
        #reform_dofset = self.ReformDofAtEachIteration
        #time_order = self.time_order

        #turb_aztec_parameters = ParameterList()
        #turb_aztec_parameters.set("AZ_solver","AZ_gmres");
        #turb_aztec_parameters.set("AZ_kspace",100);
        #turb_aztec_parameters.set("AZ_output","AZ_none");

        #turb_preconditioner_type = "ILU"
        #turb_preconditioner_parameters = ParameterList()
        #turb_overlap_level = 0
        #turb_nit_max = 1000
        #turb_linear_tol = 1e-9

        #turb_linear_solver =  AztecSolver(turb_aztec_parameters,turb_preconditioner_type,turb_preconditioner_parameters,turb_linear_tol,turb_nit_max,turb_overlap_level)
        #turb_linear_solver.SetScalingType(AztecScalingType.LeftScaling)

        #turbulence_model = TrilinosSpalartAllmarasTurbulenceModel(self.Comm,self.model_part,turb_linear_solver,self.domain_size,non_linear_tol,max_it,reform_dofset,time_order)
        #turbulence_model.AdaptForFractionalStep()
        #if(DES==True):
            #turbulence_model.ActivateDES(CDES);

        #self.solver.AddIterationStep(turbulence_model);

    ########################################################################
    def ActivateSmagorinsky(self,C):
        for elem in self.model_part.Elements:
            elem.SetValue(C_SMAGORINSKY,C)        

#################################################################################################
#################################################################################################   
def CreateSolver( model_part, config ):
    fluid_solver = IncompressibleFluidSolver( model_part, config.domain_size )
      
    ##default settings 
    fluid_solver.vel_toll = config.vel_toll
    if( hasattr(config,"vel_toll") ): fluid_solver.vel_toll = config.vel_toll
    if( hasattr(config,"press_toll") ): fluid_solver.press_toll = config.press_toll
    if( hasattr(config,"max_vel_its") ): fluid_solver.max_vel_its = config.max_vel_its
    if( hasattr(config,"max_press_its") ): fluid_solver.max_press_its = config.max_press_its
    if( hasattr(config,"time_order") ): fluid_solver.time_order = config.time_order
    if( hasattr(config,"compute_reactions") ): fluid_solver.compute_reactions = config.compute_reactions
    if( hasattr(config,"ReformDofAtEachIteration") ): fluid_solver.ReformDofAtEachIteration = config.ReformDofAtEachIteration
    if( hasattr(config,"predictor_corrector") ): fluid_solver.predictor_corrector = config.predictor_corrector
    if( hasattr(config,"echo_level") ): fluid_solver.echo_level = config.echo_level
    if( hasattr(config,"dynamic_tau") ): fluid_solver.dynamic_tau = config.dynamic_tau

    #linear solver settings
    import trilinos_linear_solver_factory
    if( hasattr(config,"pressure_linear_solver_config") ): fluid_solver.pressure_linear_solver =  trilinos_linear_solver_factory.ConstructSolver(config.pressure_linear_solver_config)
    if( hasattr(config,"velocity_linear_solver_config") ): fluid_solver.velocity_linear_solver =  trilinos_linear_solver_factory.ConstructSolver(config.velocity_linear_solver_config)
    
    #RANS or DES settings
    if( hasattr(config,"use_spalart_allmaras") ): fluid_solver.use_spalart_allmaras = config.use_spalart_allmaras
    if( hasattr(config,"use_des") ): fluid_solver.use_des = config.use_des
    if( hasattr(config,"use_spalart_allmaras")  and config.use_spalart_allmaras == True):
        if( hasattr(config,"wall_nodes") ):
            fluid_solver.wall_nodes = config.wall_nodes
        else:
            print "ATTENTION: attempting to use SpalartAllmaras without prescribig the wall position. please set the variable \"wall_nodes\" within the "
            print "config class to an appropriate list or deactivate the turbulence model"
            
    
    return fluid_solver