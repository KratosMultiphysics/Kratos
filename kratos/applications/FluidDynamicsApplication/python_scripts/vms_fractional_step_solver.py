# -*- coding: utf-8 -*-
# importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.IncompressibleFluidApplication import *
from KratosMultiphysics.FluidDynamicsApplication import *
# Check that KratosMultiphysics was imported in the main script
CheckForPreviousImport()


def AddVariables(model_part):
    model_part.AddNodalSolutionStepVariable(VELOCITY)
    model_part.AddNodalSolutionStepVariable(FRACT_VEL)
    model_part.AddNodalSolutionStepVariable(MESH_VELOCITY)
    model_part.AddNodalSolutionStepVariable(PRESSURE)
    model_part.AddNodalSolutionStepVariable(PRESSURE_OLD_IT)
    model_part.AddNodalSolutionStepVariable(PRESS_PROJ)
    model_part.AddNodalSolutionStepVariable(CONV_PROJ)
# model_part.AddNodalSolutionStepVariable(ADVPROJ)
    model_part.AddNodalSolutionStepVariable(DIVPROJ)
    model_part.AddNodalSolutionStepVariable(NODAL_AREA)
    model_part.AddNodalSolutionStepVariable(BODY_FORCE)
    model_part.AddNodalSolutionStepVariable(DENSITY)
    model_part.AddNodalSolutionStepVariable(VISCOSITY)
    model_part.AddNodalSolutionStepVariable(FLAG_VARIABLE)
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT)
    model_part.AddNodalSolutionStepVariable(IS_STRUCTURE)
    model_part.AddNodalSolutionStepVariable(REACTION)
    model_part.AddNodalSolutionStepVariable(Y_WALL)
    model_part.AddNodalSolutionStepVariable(NORMAL)
    model_part.AddNodalSolutionStepVariable(MOLECULAR_VISCOSITY)
    model_part.AddNodalSolutionStepVariable(TURBULENT_VISCOSITY)
    model_part.AddNodalSolutionStepVariable(TEMP_CONV_PROJ)
    model_part.AddNodalSolutionStepVariable(DISTANCE)
    print "variables for the vms fluid solver added correctly"


def AddDofs(model_part):

    for node in model_part.Nodes:
        # adding dofs
        node.AddDof(PRESSURE)
        node.AddDof(VELOCITY_X)
        node.AddDof(VELOCITY_Y)
        node.AddDof(VELOCITY_Z)

    print "dofs for the vms fluid solver added correctly"


class IncompressibleFluidSolver:

    def __init__(self, model_part, domain_size):

        # neighbour search
        number_of_avg_elems = 10
        number_of_avg_nodes = 10
        self.neighbour_search = FindNodalNeighboursProcess(
            model_part, number_of_avg_elems, number_of_avg_nodes)

        self.model_part = model_part
        self.domain_size = domain_size

        # assignation of parameters to be used
        self.vel_toll = 1e-20
        # 0.001
        self.press_toll = 1e-7  # 0.001;
        self.max_vel_its = 6
        self.max_press_its = 3
        self.time_order = 2
        self.CalculateReactions = False
        self.ReformDofAtEachIteration = False
        self.CalculateNormDxFlag = True
        self.laplacian_form = 1
        # 1 = laplacian, 2 = Discrete Laplacian
        self.predictor_corrector = False

        self.echo_level = 1

        # definition of the solvers
        pDiagPrecond = DiagonalPreconditioner()
#        pILUPrecond = ILU0Preconditioner()
# self.velocity_linear_solver =  BICGSTABSolver(1e-6, 5000,pDiagPrecond)
# self.pressure_linear_solver =  BICGSTABSolver(1e-9, 5000,pILUPrecond)
        self.velocity_linear_solver = ScalingSolver(
            BICGSTABSolver(1e-6, 5000, pDiagPrecond), True)
# self.pressure_linear_solver =  BICGSTABSolver(1e-3, 5000,pILUPrecond)
        # self.pressure_linear_solver =  BICGSTABSolver(1e-6,
        # 5000,pDiagPrecond)

        assume_constant_structure = True
        max_reduced_size = 1000
        self.pressure_linear_solver = ScalingSolver(
            DeflatedCGSolver(1e-6, 5000, assume_constant_structure, max_reduced_size), True)

        self.dynamic_tau = 0.001

        self.compute_reactions = False

        self.use_slip_conditions = False
        self.use_spalart_allmaras = False
        self.wall_nodes = list()
        self.use_des = False

    def Initialize(self):
        # Componentwise Builder and solver uses neighbours to set DofSet
        # (used both for PRESSURE and TURBULENT_VISCOSITY)
        (self.neighbour_search).Execute()

        self.model_part.ProcessInfo.SetValue(DYNAMIC_TAU, self.dynamic_tau)
# self.model_part.ProcessInfo.SetValue(ACTIVATE_TAU2, self.activate_tau2);

        self.domain_size = int(self.domain_size)
        self.laplacian_form = int(self.laplacian_form)
# solver_configuration =
# FractionalStepConfiguration(self.model_part,self.velocity_linear_solver,self.pressure_linear_solver,self.domain_size,self.laplacian_form
# )

        self.ReformDofAtEachIteration = bool(self.ReformDofAtEachIteration)
        self.vel_toll = float(self.vel_toll)
        self.press_toll = float(self.press_toll)
        self.max_vel_its = int(self.max_vel_its)
        self.max_press_its = int(self.max_press_its)
        self.time_order = int(self.time_order)
        self.domain_size = int(self.domain_size)
        self.predictor_corrector = bool(self.predictor_corrector)
# self.solver = FractionalStepStrategy( self.model_part,
# solver_configuration, self.ReformDofAtEachIteration, self.vel_toll,
# self.press_toll, self.max_vel_its, self.max_press_its, self.time_order,
# self.domain_size,self.predictor_corrector)

        MoveMeshFlag = False

        # check if slip conditions are defined
        if self.use_slip_conditions == False:
            for cond in self.model_part.Conditions:
                if cond.GetValue(IS_STRUCTURE) != 0.0:
                    self.use_slip_conditions = True
                    break

        # if we use slip conditions, calculate normals on the boundary
        if self.use_slip_conditions == True:
            self.normal_util = NormalCalculationUtils()
            self.normal_util.CalculateOnSimplex(
                self.model_part, self.domain_size, IS_STRUCTURE)

# self.solver = FSStrategy(self.model_part,
# self.velocity_linear_solver,
# self.pressure_linear_solver,
# MoveMeshFlag,
# self.ReformDofAtEachIteration,
# self.vel_toll,
# self.press_toll,
# self.max_vel_its,
# self.max_press_its,
# self.time_order,
# self.domain_size,
# self.predictor_corrector)
        self.solver_settings = FractionalStepSettings(self.model_part,
                                                      self.domain_size,
                                                      self.time_order,
                                                      self.use_slip_conditions,
                                                      MoveMeshFlag,
                                                      self.ReformDofAtEachIteration)
        self.solver_settings.SetEchoLevel(self.echo_level)

        self.solver_settings.SetStrategy(StrategyLabel.Velocity,
                                         self.velocity_linear_solver,
                                         self.vel_toll,
                                         self.max_vel_its)

        self.solver_settings.SetStrategy(StrategyLabel.Pressure,
                                         self.pressure_linear_solver,
                                         self.press_toll,
                                         self.max_press_its)

        if self.use_spalart_allmaras:
            for node in self.wall_nodes:
                node.SetValue(IS_VISITED, 1.0)

            distance_calculator = BodyDistanceCalculationUtils()
            if(self.domain_size == 2):
                distance_calculator.CalculateDistances2D(
                    self.model_part.Elements, DISTANCE, 100.0)
            else:
                distance_calculator.CalculateDistances3D(
                    self.model_part.Elements, DISTANCE, 100.0)

            sa_non_linear_tol = 0.001
            sa_max_it = 10

            reform_dofset = self.ReformDofAtEachIteration
            time_order = self.time_order
            pPrecond = DiagonalPreconditioner()
            turbulence_linear_solver = BICGSTABSolver(1e-9, 5000, pPrecond)

            self.solver_settings.SetTurbulenceModel(
                TurbulenceModelLabel.SpalartAllmaras,
                turbulence_linear_solver,
                sa_non_linear_tol,
                sa_max_it)

        self.solver = FSStrategy(self.model_part,
                                 self.solver_settings,
                                 self.predictor_corrector)

        self.solver.Check()

# self.solver.ApplyFractionalVelocityFixity()
        # generating the slip conditions
# self.create_slip_conditions.Execute()
# (self.solver).SetSlipProcess(self.create_slip_conditions);
# self.slip_conditions_initialized = True
# (self.solver).SetEchoLevel(self.echo_level)
        print "finished initialization of the fluid strategy"

    def Solve(self):
        if(self.ReformDofAtEachIteration == True):
# self.solver.ApplyFractionalVelocityFixity()
            (self.neighbour_search).Execute()
# self.slip_conditions_initialized = False
            if self.use_slip_conditions == True:
                self.normal_util.CalculateOnSimplex(
                    self.model_part, self.domain_size, IS_STRUCTURE)

# if(self.slip_conditions_initialized == False):
# self.create_slip_conditions.Execute()
# (self.solver).SetSlipProcess(self.create_slip_conditions);
# self.slip_conditions_initialized = True

#        print "just before solve"
#        print self.model_part.ProcessInfo.GetValue(TIME)
        (self.solver).Solve()

        if(self.compute_reactions == True):
            self.solver.CalculateReactions()  # REACTION)

# (self.create_slip_conditions).SetNormalVelocityToZero()
# (self.create_slip_conditions).ApplyEdgeConstraints()

    def Clear(self):
        (self.solver).Clear()
        self.slip_conditions_initialized = True

    def AdaptMesh(self):
        import KratosMultiphysics.MeshingApplication as KMesh
        admissible_ratio = 0.05
        max_levels = 2
        refinement_utils = KMesh.RefinementUtilities()
        if(self.domain_size == 2):
            raise "error refine in 2d not yet implemented"
        else:
            Refine = KMesh.LocalRefineTetrahedraMesh(self.model_part)
        (self.model_part).ProcessInfo[
            FRACTIONAL_STEP] = 10
        # just to be sure nothign is done
        refinement_utils.MarkForRefinement(
            ERROR_RATIO, self.model_part, admissible_ratio, max_levels)
        self.Clear()
        refine_on_reference = False
        interpolate_internal_variables = False
        Refine.LocalRefineMesh(
            refine_on_reference, interpolate_internal_variables)

        (self.neighbour_search).Execute()
        self.slip_conditions_initialized = False
        print "Refining finished"

    def WriteRestartFile(self, FileName):
        restart_file = open(FileName + ".mdpa", 'w')
        import new_restart_utilities
        new_restart_utilities.PrintProperties(restart_file)
        new_restart_utilities.PrintNodes(self.model_part.Nodes, restart_file)
        new_restart_utilities.PrintElements(
            "Fluid3D", self.model_part.Elements, restart_file)
        new_restart_utilities.PrintRestart_ScalarVariable(
            VELOCITY_X, "VELOCITY_X", self.model_part.Nodes, restart_file)
        new_restart_utilities.PrintRestart_ScalarVariable(
            VELOCITY_Y, "VELOCITY_Y", self.model_part.Nodes, restart_file)
        new_restart_utilities.PrintRestart_ScalarVariable(
            VELOCITY_Z, "VELOCITY_Z", self.model_part.Nodes, restart_file)
        new_restart_utilities.PrintRestart_ScalarVariable(
            PRESSURE, "PRESSURE", self.model_part.Nodes, restart_file)
        new_restart_utilities.PrintRestart_ScalarVariable(
            VISCOSITY, "VISCOSITY", self.model_part.Nodes, restart_file)
        new_restart_utilities.PrintRestart_ScalarVariable(
            DENSITY, "DENSITY", self.model_part.Nodes, restart_file)
        restart_file.close()

    def ActivateSmagorinsky(self, C):
        for elem in self.model_part.Elements:
            elem.SetValue(C_SMAGORINSKY, C)

    def ActivateSpalartAllmaras(self, wall_nodes, DES, CDES=1.0):
        self.use_spalart_allmaras = True
        self.wall_nodes = wall_nodes
        # import KratosMultiphysics.FluidDynamicsApplication as KCFD
        # for node in wall_nodes:
            # node.SetValue(IS_VISITED,1.0)

        # distance_calculator = BodyDistanceCalculationUtils()
        # distance_calculator.CalculateDistances2D(self.model_part.Elements,DISTANCE,100.0)

        # non_linear_tol = 0.001
        # max_it = 10
        # reform_dofset = self.ReformDofAtEachIteration
        # time_order = self.time_order
        # pPrecond = DiagonalPreconditioner()
        # turbulence_linear_solver =  BICGSTABSolver(1e-20, 5000,pPrecond)
        # turbulence_model = KCFD.SpalartAllmarasTurbulenceModel(self.model_part,turbulence_linear_solver,self.domain_size,non_linear_tol,max_it,reform_dofset,time_order);
# turbulence_model.AdaptForFractionalStep()
        # if(DES==True):
            # turbulence_model.ActivateDES(CDES);

        # self.solver.AddIterationStep(turbulence_model);
        
        
        
        
        
        
        
       
#################################################################################################
#################################################################################################
def AddDofUtility( model_part, config ):
    AddDofs( model_part )
    
def AddVariablesUtility( model_part, config ):
    AddVariables( model_part )
    
def CreateSolver( model_part, config ):
    fluid_solver = IncompressibleFluidSolver( model_part, config.domain_size )
    
    #import config_reader
    #config_reader.AssignRequiredValue( fluid_solver.vel_toll , config, "vel_toll")
    #config_reader.AssignRequiredValue( fluid_solver.press_toll , config, "press_toll")
    #config_reader.AssignRequiredValue( fluid_solver.max_vel_its , config, "max_vel_its")
    #config_reader.AssignRequiredValue( fluid_solver.max_press_its , config, "max_press_its")
    #config_reader.AssignRequiredValue( fluid_solver.time_order , config, "time_order")
    #config_reader.AssignRequiredValue( fluid_solver.compute_reactions , config, "compute_reactions")
    #config_reader.AssignRequiredValue( fluid_solver.predictor_corrector , config, "predictor_corrector")
    
    #config_reader.AssignOptionalValue( fluid_solver.ReformDofAtEachIteration , config, "ReformDofAtEachIteration")
    #config_reader.AssignOptionalValue( fluid_solver.echo_level , config, "echo_level")
    #config_reader.AssignOptionalValue( fluid_solver.dynamic_tau , config, "dynamic_tau")
    
    
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
    import linear_solver_factory
    if( hasattr(config,"pressure_linear_solver_config") ): fluid_solver.pressure_linear_solver =  linear_solver_factory.ConstructSolver(config.pressure_linear_solver_config)
    if( hasattr(config,"velocity_linear_solver_config") ): fluid_solver.velocity_linear_solver =  linear_solver_factory.ConstructSolver(config.velocity_linear_solver_config)

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
    