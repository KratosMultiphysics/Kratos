from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.IncompressibleFluidApplication import *
from KratosMultiphysics.ThermoMechanicalApplication import *
from KratosMultiphysics.MeshingApplication import *
from KratosMultiphysics.FluidDynamicsApplication import *
from KratosMultiphysics.ExternalSolversApplication import *
# Check that KratosMultiphysics was imported in the main script
CheckForPreviousImport()

import levelset_solver

# settings for the convection solver
distance_settings = ConvectionDiffusionSettings()
distance_settings.SetUnknownVariable(DISTANCE)
distance_settings.SetConvectionVariable(VELOCITY)
distance_settings.SetMeshVelocityVariable(MESH_VELOCITY)
# distance_settings.SetVolumeSourceVariable(HEAT_FLUX)
# distance_settings.SetDiffusionVariable(ARRHENIUSAUX)
# For level set solver rho and C are assigned to 1
distance_settings.SetDensityVariable(DISTANCE)


def AddVariables(model_part):
    model_part.AddNodalSolutionStepVariable(VELOCITY)
    model_part.AddNodalSolutionStepVariable(ACCELERATION)
    model_part.AddNodalSolutionStepVariable(MESH_VELOCITY)
    model_part.AddNodalSolutionStepVariable(PRESSURE)
    model_part.AddNodalSolutionStepVariable(IS_FLUID)
    model_part.AddNodalSolutionStepVariable(IS_STRUCTURE)
    model_part.AddNodalSolutionStepVariable(IS_FREE_SURFACE)
    model_part.AddNodalSolutionStepVariable(IS_AIR_EXIT)
    model_part.AddNodalSolutionStepVariable(IS_INTERFACE)
    model_part.AddNodalSolutionStepVariable(IS_BOUNDARY)
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT)
    model_part.AddNodalSolutionStepVariable(VISCOSITY)
    model_part.AddNodalSolutionStepVariable(DENSITY)
    model_part.AddNodalSolutionStepVariable(BODY_FORCE)
    model_part.AddNodalSolutionStepVariable(NODAL_AREA)
    model_part.AddNodalSolutionStepVariable(NODAL_H)
    model_part.AddNodalSolutionStepVariable(THAWONE)
    model_part.AddNodalSolutionStepVariable(THAWTWO)
    model_part.AddNodalSolutionStepVariable(FLAG_VARIABLE)
    model_part.AddNodalSolutionStepVariable(NORMAL)
    model_part.AddNodalSolutionStepVariable(DISTANCE)
    model_part.AddNodalSolutionStepVariable(IS_SLIP)
    model_part.AddNodalSolutionStepVariable(PRESSURES)
    model_part.AddNodalSolutionStepVariable(VELOCITIES)
    model_part.AddNodalSolutionStepVariable(MATERIAL)
    model_part.AddNodalSolutionStepVariable(LAST_AIR)
    model_part.AddNodalSolutionStepVariable(NODAL_MASS)
    model_part.AddNodalSolutionStepVariable(NODAL_PAUX)
    model_part.AddNodalSolutionStepVariable(Y_WALL)
    model_part.AddNodalSolutionStepVariable(FILLTIME)
    model_part.AddNodalSolutionStepVariable(MAX_VEL)
    model_part.AddNodalSolutionStepVariable(WET_VOLUME) #
    model_part.AddNodalSolutionStepVariable(MOULD_INNER_TEMPERATURE)
    # variables needed for the distance solver
    levelset_solver.AddVariables(model_part, distance_settings)

    print("variables for the MONOLITHIC_SOLVER_EULERIAN added correctly")


def AddDofs(model_part):
    for node in model_part.Nodes:
        # adding dofs
        node.AddDof(VELOCITY_X)
        node.AddDof(VELOCITY_Y)
        node.AddDof(VELOCITY_Z)
        node.AddDof(PRESSURE)

        #node.AddDof(MESH_VELOCITY_X)
        #node.AddDof(MESH_VELOCITY_Y)
        #node.AddDof(MESH_VELOCITY_Z)

    levelset_solver.AddDofs(model_part, distance_settings)
    print("dofs for the monolithic solver added correctly")


class MonolithicSolver:
    #

    def __init__(self, model_part, domain_size,linear_solver_iterations=300, linear_solver_tolerance=1e-5):

        self.model_part = model_part
        self.domain_size = domain_size

        self.alpha = -0.0
        self.move_mesh_strategy = 0
        # self.time_scheme = ResidualBasedPredictorCorrectorVelocityBossakSchemeDPGEnriched(
            # self.alpha, self.move_mesh_strategy, self.domain_size)

        self.time_scheme = ResidualBasedPredictorCorrectorBDFSchemeTurbulent(self.domain_size)
        self.time_scheme.Check(self.model_part)

        # self.time_scheme = ResidualBasedPredictorCorrectorVelocityBossakScheme( self.alpha,self.move_mesh_strategy )
        # definition of the solvers
        # self.linear_solver =  SkylineLUFactorizationSolver()
# self.linear_solver =SuperLUSolver()
# self.linear_solver = MKLPardisoSolver()
        # pPrecond = DiagonalPreconditioner()
# pPrecond = ILU0Preconditioner()
         # self.linear_solver =  BICGSTABSolver(1e-6, 5000,pPrecond)
        # gmres_size = 30
        # ilu_level_of_fill = 2
        # tol = 1e-5
        # verbosity = 0
        # self.linear_solver = PastixSolver(tol,gmres_size,ilu_level_of_fill,verbosity,False)
        # self.linear_solver = PastixSolver(verbosity,False)
        # new solvers
        self.gmres_size = 200
        self.iterations = linear_solver_iterations #400 # Ojo, antes 200
        self.tol = linear_solver_tolerance #1e-5 #Before 1e-5
        self.verbosity = 1
        self.linear_solver = ScalingSolver(AMGCLSolver(
            AMGCLSmoother.ILU0,
            AMGCLIterativeSolverType.BICGSTAB_WITH_GMRES_FALLBACK, #AMGCLIterativeSolverType.GMRES, #BICGSTAB_WITH_GMRES_FALLBACK,
            self.tol,
            self.iterations,
            self.verbosity,
            self.gmres_size),True)
        print("#####################################")
        print("#### LINEAR SOLVER               ####")
        print("#####################################")
        print(self.linear_solver)
        print("#####################################")
        print("A VER ")
#### PRUEBA a ver 

        # definition of the convergence criteria
        self.rel_vel_tol = 1e-5
        self.abs_vel_tol = 1e-7
        self.rel_pres_tol = 1e-5
        self.abs_pres_tol = 1e-7

        self.dynamic_tau_levelset = 0.01
        self.dynamic_tau_fluid = 0.01
        self.oss_switch = 0

        # non newtonian setting
        self.regularization_coef = 1000

        self.max_iter = 30

        # default settings
        self.echo_level = 0
        self.CalculateReactionFlag = True
        self.ReformDofSetAtEachStep = True
        self.CalculateNormDxFlag = True
        self.MoveMeshFlag = False
        self.volume_correction = True
        self.vol_cr_step = 5

# print "Construction monolithic solver finished"

        # Creat Lavel_set solver
        # construct the model part
        if(domain_size == 2):
            raise "error, still not implemented in 2D"
            conv_elem = "SUPGConv2D"
            conv_cond = "Condition2D"
        else:
            conv_elem = "SUPGConv3D"
            conv_cond = "Condition3D"
        self.level_set_model_part = ModelPart("level_set_model_part")
        self.conv_generator = ConnectivityPreserveModeler()
        (self.conv_generator).GenerateModelPart(self.model_part,
                                                self.level_set_model_part, conv_elem, conv_cond)
        #(ParallelFillCommunicator(self.level_set_model_part)).Execute();

        # constructing the convection solver for the distance
        self.level_set_solver = levelset_solver.Solver(self.level_set_model_part,domain_size,distance_settings)
        self.level_set_solver.max_iter = 8

        #
        # properties of the two fluids
        self.rho1 = 2400.0  # applied on the negative part of the domain 1000.0
        self.conductivity1 = 1.0

        self.rho2 = 1.0  # applied to the positive part of the domain#1.0
        self.conductivity2 = 1.0

        self.mu = 3.0e-3
        self.divergence_clearance_performed = False
        #

        # Distance utilities

         #
        if(self.domain_size == 2):
            self.redistance_utils = ParallelDistanceCalculator2D()
        else:
            self.redistance_utils = ParallelDistanceCalculator3D()

        self.redistance_frequency = 1
        self.max_edge_size = self.redistance_utils.FindMaximumEdgeSize(
            self.level_set_model_part)
        self.max_distance = self.max_edge_size * 5.0
        self.max_levels = 25  # self.max_distance/self.min_edge_size

        self.max_ns_iterations = 8
        self.internal_step_counter = 1

        # Slip condition
        self.use_slip_conditions = False

        # volume correction
        self.volume_correction_switch = True
        self.negative_volume_correction=True

        # element size
        self.maxmin = []
        ParticleLevelSetUtils3D().FindMaxMinEdgeSize(
            self.level_set_model_part, self.maxmin)

    #

    def ApplyFluidProperties(self):
        # apply density
        mu1 = 1.0 * self.mu / self.rho1
        # mu1 = self.mu
        # mu2 = 0.01*self.mu/self.rho2
        mu2 = mu1
        BiphasicFillingUtilities().ApplyFluidProperties(self.model_part, mu1, self.rho1, mu2, self.rho2)
# for node in self.model_part.Nodes:
# dist = node.GetSolutionStepValue(DISTANCE)
# if(dist < 0):
# node.SetSolutionStepValue(DENSITY,0,self.rho1)
# node.SetSolutionStepValue(VISCOSITY,0,mu1)
# else:
# node.SetSolutionStepValue(DENSITY,0,self.rho2)
# node.SetSolutionStepValue(VISCOSITY,0,mu2)
    #

    def Initialize(self):
        # creating the solution strategy
        self.conv_criteria = VelPrCriteria(self.rel_vel_tol, self.abs_vel_tol,
                                           self.rel_pres_tol, self.abs_pres_tol)
        builder_and_solver = ResidualBasedBlockBuilderAndSolver(
            self.linear_solver)
        # self.solver = ResidualBasedNewtonRaphsonStrategy(self.model_part,self.time_scheme,self.linear_solver,self.conv_criteria,builder_and_solver,self.max_iter,self.CalculateReactionFlag, self.ReformDofSetAtEachStep,self.MoveMeshFlag)

        self.solver = ResidualBasedNewtonRaphsonStrategy(
            self.model_part,
            self.time_scheme,
            self.linear_solver,
            self.conv_criteria,
            builder_and_solver,
            self.max_iter,
            self.CalculateReactionFlag,
            self.ReformDofSetAtEachStep,
            self.MoveMeshFlag)
        (self.solver).SetEchoLevel(self.echo_level)
        #print(">>>>>>>>>>>>>>> OSS_SWITCH = ", self.oss_switch)
        self.model_part.ProcessInfo.SetValue(
            DYNAMIC_TAU, self.dynamic_tau_fluid)
        self.model_part.ProcessInfo.SetValue(OSS_SWITCH, self.oss_switch)

        # LEvel_set solver initialization
        self.level_set_solver.dynamic_tau = self.dynamic_tau_levelset
        
        
        self.redistance_utils.CalculateDistances(self.model_part,DISTANCE, NODAL_AREA,self.max_levels,self.max_distance)
        
        # self.redistance_utils.CalculateInterfacePreservingDistances(self.model_part,DISTANCE,NODAL_AREA,self.max_levels,self.max_distance)
        self.level_set_solver.linear_solver = AMGCLSolver(
            AMGCLSmoother.ILU0,
            AMGCLIterativeSolverType.GMRES,
            1e-6, #self.tol,
            200,
            self.verbosity,
            self.gmres_size)
        self.level_set_solver.Initialize()

        self.ApplyFluidProperties()

        self.next_redistance = self.redistance_frequency
        #
        # FOR SLIP
        #
        # Manullay assign!
        for cond in self.model_part.Conditions:
            cond.SetValue(IS_STRUCTURE, 1.0)
        # if we use slip conditions, calculate normals on the boundary
        if (self.use_slip_conditions):
            (FindConditionsNeighboursProcess(
                self.model_part, 3, 20)).ClearNeighbours()
            (FindConditionsNeighboursProcess(self.model_part, 3, 20)).Execute()
            self.normal_util = NormalCalculationUtils()
            self.normal_util.CalculateOnSimplex(
                self.model_part,
                self.domain_size,
                IS_STRUCTURE, 0, 35.0)  # ,0.0),180) #35.0)  # ,0.0,35.0

        # saving inlet nodes
        self.inlet_nodes = []
        for cond in self.model_part.Conditions:
            if(cond.GetValue(IS_INLET) > 0):
                for node in cond.GetNodes():
                    self.inlet_nodes.append(node)
                    
                    
                    
        #import velocity_convection_utility
        #self.velocity_prediction = velocity_convection_utility.VelocityConvectionUtility(self.model_part)

    #
    def DoRedistance(self):
        # redistance if required
        self.redistance_utils.CalculateDistances(self.model_part,DISTANCE,NODAL_AREA, self.max_levels, self.max_distance)
        # self.redistance_utils.CalculateInterfacePreservingDistances(self.model_part,DISTANCE,NODAL_AREA,self.max_levels,self.max_distance)

     #
    def ConvectDistance(self):
        #self.level_set_model_part.ProcessInfo = self.model_part.ProcessInfo
        (self.level_set_model_part.ProcessInfo).SetValue(CONVECTION_DIFFUSION_SETTINGS, distance_settings)
        (self.level_set_model_part.ProcessInfo).SetValue(DYNAMIC_TAU, self.dynamic_tau_levelset)  # self.dynamic_tau
        (self.level_set_solver).Solve()
        BiphasicFillingUtilities().DistanceFarRegionCorrection(self.model_part,self.max_distance)
     #
      #

    def Solve(self, step):
        # at the beginning of the calculations do a div clearance step
        if(self.divergence_clearance_performed == False):
            for node in self.model_part.Nodes:
                node.SetSolutionStepValue(DISTANCE,1,node.GetSolutionStepValue(DISTANCE))
            self.divergence_clearance_performed = True

        if(step > 3):
            (self.solver).Predict()

        Timer.Start("ConvectDistance")
        # convect distance function

        self.ConvectDistance()

        # recompute distance function as needed
        Timer.Start("DoRedistance")
        if(self.internal_step_counter >= self.next_redistance):
            #net_volume = self.model_part.ProcessInfo[NET_INPUT_MATERIAL]
            #BiphasicFillingUtilities().VolumeCorrection(self.model_part, net_volume, self.max_edge_size)
            
            #ensure that inlet nodes are still wet
            for node in self.inlet_nodes:
                if( node.GetSolutionStepValue(DISTANCE) > 0.0):
                    node.SetSolutionStepValue(DISTANCE,0,  -0.01*self.max_edge_size)
            self.DoRedistance()
        Timer.Stop("DoRedistance")

        if(self.volume_correction_switch and step > self.vol_cr_step):
            net_volume = self.model_part.ProcessInfo[NET_INPUT_MATERIAL]
            BiphasicFillingUtilities().VolumeCorrection(self.model_part, net_volume, self.max_edge_size,self.negative_volume_correction)
        Timer.Start("ApplyFluidProperties")
        self.ApplyFluidProperties()
        BiphasicFillingUtilities().ViscosityBasedSolidification(self.model_part,100.0)
        #self.IncreaseCSmagToSOlidify(50.0)
        Timer.Stop("ApplyFluidProperties")
        # Recompute normals if necessary
# if(self.ReformDofSetAtEachStep == True):
# if self.use_slip_conditions == True:
# self.normal_util.CalculateOnSimplex(self.model_part,self.domain_size,IS_STRUCTURE,0.0,35.0)#,0.0,35.0

        Timer.Start("self.solve")

        if(step > 3):
            (self.solver).Predict()
            
            #self.velocity_prediction.PredictVelocity()
            
           
           
        #ActivationUtilities().ActivateElementsAndConditions( self.model_part, DISTANCE, self.max_distance, True) 
        (self.solver).Solve()
        self.internal_step_counter += 1
        Timer.Stop("self.solve")
    #

    def SetEchoLevel(self, level):
        (self.solver).SetEchoLevel(level)

    #
##    def IncreaseCSmagToSOlidify(self,val):
##        for elem in self.model_part.Elements:
##            max_alpha = 0.0
##            for node in elem.GetNodes():
##                alpha = node.GetSolutionStepValue(DP_ALPHA1)
##                if(alpha > max_alpha):
##                    max_alpha = alpha
##            if(max_alpha > 0.0):
##                new_csmag = 0.45 + max_alpha*val #0.45 is default CSmag of the calculation
##                elem.SetValue(C_SMAGORINSKY, new_csmag )
