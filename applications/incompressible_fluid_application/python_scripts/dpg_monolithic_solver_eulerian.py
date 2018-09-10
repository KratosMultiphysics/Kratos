from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.IncompressibleFluidApplication import *
from KratosMultiphysics.ThermoMechanicalApplication import *
#from KratosMultiphysics.MeshingApplication import *
# import KratosMultiphysics.MeshingApplication as Meshing
from KratosMultiphysics.FluidDynamicsApplication import *
#from KratosMultiphysics.ExternalSolversApplication import *
import KratosMultiphysics.ExternalSolversApplication as ExternalSolver

from KratosMultiphysics.Click2CastApplication import *
import linear_solver_factory

import sys
import traceback
import math
import time as time_measure
import os
class ComputeFlowLength1:
    """This Calss is used to compute a very raw stimation of the flow length using that the variation in distance is equal to the distance of  """
    def __init__(self,ModelPart):
        self.ModelPart=ModelPart
        # Set FlowLength1 to 0
        for node in self.ModelPart.Nodes:
            node.SetValue(FLOW_LENGTH,0.0)
    def Execute(self):
        for node in self.ModelPart.Nodes:
            d1=node.GetSolutionStepValue(DISTANCE,1)
            if(d1<0.0):
                d0=node.GetSolutionStepValue(DISTANCE)
                # Read the original Value
                flength_old=node.GetValue(FLOW_LENGTH)
                delta_distance=-(d0-d1)
                node.SetValue(FLOW_LENGTH,flength_old+delta_distance)



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
    #model_part.AddNodalSolutionStepVariable(IS_FREE_SURFACE)
    model_part.AddNodalSolutionStepVariable(IS_AIR_EXIT)
    model_part.AddNodalSolutionStepVariable(IS_INTERFACE)
    model_part.AddNodalSolutionStepVariable(IS_BOUNDARY)
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT)
    model_part.AddNodalSolutionStepVariable(VISCOSITY)
    model_part.AddNodalSolutionStepVariable(DENSITY)
    model_part.AddNodalSolutionStepVariable(BODY_FORCE)
    model_part.AddNodalSolutionStepVariable(NODAL_AREA)
    #model_part.AddNodalSolutionStepVariable(NODAL_H)
    model_part.AddNodalSolutionStepVariable(THAWONE)
    model_part.AddNodalSolutionStepVariable(THAWTWO)
    model_part.AddNodalSolutionStepVariable(FLAG_VARIABLE)
    model_part.AddNodalSolutionStepVariable(NORMAL)
    model_part.AddNodalSolutionStepVariable(DISTANCE)
    model_part.AddNodalSolutionStepVariable(IS_SLIP)
    model_part.AddNodalSolutionStepVariable(PRESSURES)
    #model_part.AddNodalSolutionStepVariable(MATERIAL)
    #model_part.AddNodalSolutionStepVariable(LAST_AIR)
    #model_part.AddNodalSolutionStepVariable(NODAL_MASS)
    model_part.AddNodalSolutionStepVariable(NODAL_PAUX)
    model_part.AddNodalSolutionStepVariable(Y_WALL)
    #model_part.AddNodalSolutionStepVariable(FILLTIME)
    #model_part.AddNodalSolutionStepVariable(MAX_VEL)
    #model_part.AddNodalSolutionStepVariable(WET_VOLUME) #
    #model_part.AddNodalSolutionStepVariable(DISTANCE) #

    print("variables for the MONOLITHIC_SOLVER_EULERIAN added correctly")


def AddDofs(model_part):
    for node in model_part.Nodes:
        # adding dofs
        node.AddDof(VELOCITY_X)
        node.AddDof(VELOCITY_Y)
        node.AddDof(VELOCITY_Z)
        node.AddDof(PRESSURE)
    print("dofs for the monolithic solver added correctly")


class MonolithicSolver:

    def __init__(self, model_part, domain_size,fluid_linear_solver_settings,redistance_settings,ConvectionType):
        Timer.Start("DPG_Monolithic_constructor")
        if(ConvectionType=="Classic"):
            self.use_distance_convection_process=False
        else:
            self.use_distance_convection_process=True
        # DISTANCE CONVECTION SELECTION
        if(self.use_distance_convection_process==False):
            import levelset_solver
            levelset_solver.AddVariables(model_part, distance_settings)
            levelset_solver.AddDofs(model_part, distance_settings)
        else:
            pass
            #model_part.AddNodalSolutionStepVariable(DISTANCE)
        #self.compute_flow_length_object=ComputeFlowLength1(model_part)
    #def __init__(self, model_part, domain_size,linear_solver_iterations=300, linear_solver_tolerance=1e-5,dynamic_tau_levelset=0.01):
        self.fluid_linear_solver_settings=fluid_linear_solver_settings
        self.redistance_settings=redistance_settings
        self.redistance_type=redistance_settings.redistance_type
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
        self.gmres_size = fluid_linear_solver_settings.gmres_size
        self.iterations = fluid_linear_solver_settings.linear_solver_iterations #400 # Ojo, antes 200
        self.tol = fluid_linear_solver_settings.linear_solver_tolerance #1e-5 #Before 1e-5
        self.verbosity = fluid_linear_solver_settings.verbosity
        if(fluid_linear_solver_settings.linear_solver_type=="BICGSTAB2"):
            #self.linear_solver = ExternalSolver.ScalingSolver(AMGCLSolver(
            self.linear_solver = ScalingSolver(ExternalSolver.AMGCLSolver(
                ExternalSolver.AMGCLSmoother.ILU0,
                ExternalSolver.AMGCLIterativeSolverType.BICGSTAB2, #AMGCLIterativeSolverType.GMRES, #BICGSTAB_WITH_GMRES_FALLBACK,
                fluid_linear_solver_settings.linear_solver_tolerance,
                fluid_linear_solver_settings.linear_solver_iterations,
                fluid_linear_solver_settings.verbosity,
                fluid_linear_solver_settings.gmres_size),True)
        else:
            self.linear_solver = ScalingSolver(ExternalSolver.AMGCLSolver(
                ExternalSolver.AMGCLSmoother.ILU0,
                ExternalSolver.AMGCLIterativeSolverType.BICGSTAB_WITH_GMRES_FALLBACK, #AMGCLIterativeSolverType.GMRES, #BICGSTAB_WITH_GMRES_FALLBACK,
                fluid_linear_solver_settings.linear_solver_tolerance,
                fluid_linear_solver_settings.linear_solver_iterations,
                fluid_linear_solver_settings.verbosity,
                fluid_linear_solver_settings.gmres_size),True)
        # definition of the convergence criteria
        self.rel_vel_tol = 1e-5
        self.abs_vel_tol = 1e-7
        self.rel_pres_tol = 1e-5
        self.abs_pres_tol = 1e-7

        self.dynamic_tau_levelset =  fluid_linear_solver_settings.dynamic_tau_levelset
        self.dynamic_tau_fluid = fluid_linear_solver_settings.dynamic_tau_fluid #0.01
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
        # Data for constructing the mass loss
        self.AcumulatedMassLossInWalls=0.0
        self.AcumulatedMassLossInWetWalls=0.0
        self.AcumulatedConvectedMass=0.0
        self.AcumulatedLossInRedistance=0.0

# print "Construction monolithic solver finished"

        # Creat Lavel_set solver
        # construct the model part
        if(self.use_distance_convection_process==False):
            if(domain_size == 2):
                raise "error, still not implemented in 2D"
                conv_elem = "SUPGConv2D"
                conv_cond = "Condition2D"
            else:
                conv_elem = "SUPGConv3D" #"SUPGConvLevelSet"#"SUPGConv3D"
                conv_cond = "Condition3D"
            self.level_set_model_part = ModelPart("level_set_model_part")
            self.conv_generator = ConnectivityPreserveModeler()
            (self.conv_generator).GenerateModelPart(self.model_part,self.level_set_model_part, conv_elem, conv_cond)
            # constructing the convection solver for the distance
            self.level_set_solver = levelset_solver.Solver(self.level_set_model_part,domain_size,distance_settings)
            self.level_set_solver.max_iter = 8
        else: #### The Distance Convection will be used
            max_cfl = 5.0;
            ## For using linear solver GMRES #AMGCLIterativeSolverType.GMRES, #BICGSTAB_WITH_GMRES_FALLBACK,
            #linear_solver = ScalingSolver(AMGCLSolver( AMGCLSmoother.ILU0, AMGCLIterativeSolverType.GMRES, 1e-9, 50, fluid_linear_solver_settings.verbosity, 50),True)
            self.convection_tolerance=1e-9
            self.convection_max_number_iterations=100
            self.convection_gmres_size=25
            linear_solver = ScalingSolver(ExternalSolver.AMGCLSolver( ExternalSolver.AMGCLSmoother.ILU0, ExternalSolver.AMGCLIterativeSolverType.BICGSTAB_WITH_GMRES_FALLBACK,self.convection_tolerance, self.convection_max_number_iterations, fluid_linear_solver_settings.verbosity, self.convection_gmres_size),True)

            self.levelset_convection_process = LevelSetConvectionProcess3D(DISTANCE,self.model_part, linear_solver, max_cfl)
     #
        # properties of the two fluids
        self.rho1 = 2400.0  # applied on the negative part of the domain 1000.0
        self.conductivity1 = 1.0

        self.rho2 = 1.0  # applied to the positive part of the domain#1.0
        self.conductivity2 = 1.0

        self.mu = 3.0e-3
        self.divergence_clearance_performed = False

        ##########################################
        ##### Compute Max Edge Size       ########
        ##########################################

        if(self.domain_size == 2):
            self.max_edge_size = ParallelDistanceCalculator2D().FindMaximumEdgeSize(self.model_part)
        else:
            self.max_edge_size = ParallelDistanceCalculator3D().FindMaximumEdgeSize(self.model_part)

        self.average_edge_size = BiphasicFillingUtilities().ComputePartAvgh(self.model_part)
        # self.max_distance = self.max_edge_size * self.redistance_settings.max_distance_factor #Ojo antes 5.0

        ##########################################
        ##### Initialize Counters         ########
        ##########################################
        self.internal_step_counter = 1
        self.redistance_frequency = self.redistance_settings.redistance_frequency

        ##########################################
        ##### OLD REDISTANCE              ########
        ##########################################


        # Distance utilities
        if(self.redistance_type=="Old"):
            print("performing Old Redistance")
            if(self.domain_size == 2):
                self.distance_calculator = ParallelDistanceCalculator2D()
            else:
                self.distance_calculator = ParallelDistanceCalculator3D()

            self.max_edge_size = self.distance_calculator.FindMaximumEdgeSize(self.model_part)
            self.max_distance = self.max_edge_size * self.redistance_settings.max_distance_factor #Ojo antes 5.0
            self.max_levels = self.redistance_settings.max_levels #Ojo antes 25

            self.max_ns_iterations =self.redistance_settings.max_ns_iterations #Ojo antes 8
        else:
            print("performing New Redistance")
            for cond in self.model_part.Conditions:
                for node in cond.GetNodes():
                    node.Set(BOUNDARY,True)
            distance_calculator_aux = ParallelDistanceCalculator3D()
            self.max_edge_size = distance_calculator_aux.FindMaximumEdgeSize(self.model_part)
            distance_linear_solver=linear_solver_factory.ConstructSolver(self.redistance_settings)
            self.distance_calculator=VariationalDistanceCalculationProcess3D(self.model_part,distance_linear_solver,self.redistance_settings.redistance_iterations)
            self.max_distance = self.max_edge_size * self.redistance_settings.max_distance_factor
        # Slip condition
        self.use_slip_conditions = False

        # volume correction
        self.volume_correction_switch = True
        self.negative_volume_correction=False
        self.volume_corrector_type="Old"

        # element size
        self.maxmin = []
        ParticleLevelSetUtils3D().FindMaxMinEdgeSize(self.model_part, self.maxmin)
        # Variables needed for computing the Efficiency of the Injection
        self.OldNetInletVolume=0.0
        self.OldWetVolume=0.0

        self.distance_utilities=DistanceUtilities()
        self.volume_corrector=CorrectVolume(self.model_part,self.average_edge_size ,self.max_edge_size)
        Timer.Stop("DPG_Monolithic_constructor")
    #

    def ApplyFluidProperties(self):
        # apply density
        mu1 = 1.0 * self.mu / self.rho1
        # mu1 = self.mu
        # mu2 = 0.01*self.mu/self.rho2
        mu2 = mu1
        BiphasicFillingUtilities().ApplyFluidProperties(self.model_part, mu1, self.rho1, mu2, self.rho2)

    def Initialize(self):
        Timer.Start("DPG_Monolithic_initialize")
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

        self.model_part.ProcessInfo.SetValue(
            DYNAMIC_TAU, self.dynamic_tau_fluid)
        self.model_part.ProcessInfo.SetValue(OSS_SWITCH, self.oss_switch)


        #Calculate Distances
        self.DoRedistance()

        # LEvel_set solver initialization
        if(self.use_distance_convection_process==False):
            self.level_set_solver.dynamic_tau = self.dynamic_tau_levelset
            self.level_set_solver.linear_solver = ExternalSolver.AMGCLSolver(
                ExternalSolver.AMGCLSmoother.ILU0,
                ExternalSolver.AMGCLIterativeSolverType.GMRES,
                1e-6, #self.tol,
                200,
                self.verbosity,
                self.gmres_size)
            self.level_set_solver.Initialize()
        # End of Levelset SOlver Definition
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
            (FindConditionsNeighboursProcess(self.model_part, 3, 20)).ClearNeighbours()
            (FindConditionsNeighboursProcess(self.model_part, 3, 20)).Execute()
            self.normal_util = NormalCalculationUtils()
            ## If needed we swap normals - Mesher is Getting Normals Inside
            self.normal_util.CalculateOnSimplex(self.model_part,self.domain_size, IS_STRUCTURE, 0, 35.0)  # ,0.0),180) #35.0)  # ,0.0,35.0
            #self.normal_util.CalculateOnSimplex(self.model_part,self.domain_size)  # ,0.0),180) #35.0)  # ,0.0,35.0
        # saving inlet nodes
        self.inlet_nodes = []
            # Memory profiler
        for cond in self.model_part.Conditions:
            if(cond.GetValue(IS_INLET) > 0):
                for node in cond.GetNodes():
                    self.inlet_nodes.append(node)
        #import velocity_convection_utility
        #self.velocity_prediction = velocity_convection_utility.VelocityConvectionUtility(self.model_part)
        Timer.Stop("DPG_Monolithic_initialize")
    #
    def DoRedistance(self):
        # redistance if required
        print("")
        print("* Recomputing Free Surface interface")
        if(self.redistance_type=="Old"):
            print("    Performing Old Redistance")
            self.distance_calculator.CalculateDistances(self.model_part,DISTANCE,NODAL_AREA, self.max_levels, self.max_distance)
            #self.distance_calculator.CalculateInterfacePreservingDistances(self.model_part,DISTANCE,NODAL_AREA, self.max_levels, self.max_distance)
        else:
            print("    Performing New Redistance")
            self.distance_calculator.Execute()

    def ConvectDistance(self):
        print("")
        print("* Convecting Free Surface")
        if(self.use_distance_convection_process==False): # Classical Way
            print( "    Convecting Distance in the Classical Way")
            (self.level_set_model_part.ProcessInfo).SetValue(CONVECTION_DIFFUSION_SETTINGS, distance_settings)
            (self.level_set_model_part.ProcessInfo).SetValue(DYNAMIC_TAU, self.dynamic_tau_levelset)  # self.dynamic_tau
            (self.level_set_solver).Solve()
        else: #New Process way
            print( "    Convecting Distance in the New Way")
            self.levelset_convection_process.Execute()
        BiphasicFillingUtilities().DistanceFarRegionCorrection(self.model_part,self.max_distance)
    def CorrectVolume(self):
        BiphasicFillingUtilities().ComputeNetInletVolume(self.model_part)
        if(self.volume_correction_switch and self.step > self.vol_cr_step):
            net_volume = self.model_part.ProcessInfo[NET_INPUT_MATERIAL]
            print("")
            print("* Performing Volume Correction")
            if(self.volume_corrector_type=="New"):
                self.volume_corrector.Execute(net_volume)
            else:
                BiphasicFillingUtilities().VolumeCorrection(self.model_part, net_volume, self.average_edge_size,self.negative_volume_correction)
        else:
            self.model_part.ProcessInfo[WET_VOLUME] = BiphasicFillingUtilities().ComputeWetVolume(self.model_part) #fluid_model_part.ProcessInfo[WET_VOLUME]
    def Solve(self, step):
        Timer.Start("DPG_Monolithic_Solve")
        print("")
        print("******************************")
        print("*     Solving for Fluid      *")
        print("******************************")
        self.step=step
        # at the beginning of the calculations do a div clearance step
        if(self.divergence_clearance_performed == False):
            for node in self.model_part.Nodes:
                node.SetSolutionStepValue(DISTANCE,1,node.GetSolutionStepValue(DISTANCE))
            self.divergence_clearance_performed = True

        # Now we store the Old Values
        self.OldNetInletVolume=self.model_part.ProcessInfo[NET_INPUT_MATERIAL]
        self.OldWetVolume=self.model_part.ProcessInfo[WET_VOLUME]
        #Timer.Start("DPG_Monolithic_Solve_ComputeWetVolume")
        #WetVolumeBeforeConvecting=BiphasicFillingUtilities().ComputeWetVolume(self.model_part)
        #Timer.Stop("DPG_Monolithic_Solve_ComputeWetVolume")
        Timer.Start("DPG_Monolithic_Solve_Predict")
        if(step > 3):
            (self.solver).Predict()
        Timer.Stop("DPG_Monolithic_Solve_Predict")
        Timer.Start("DPG_Monolithic_Solve_ConvectDistance")
        # convect distance function
        self.ConvectDistance()
        Timer.Stop("DPG_Monolithic_Solve_ConvectDistance")
        # recompute distance function as needed

        #Timer.Start("DPG_Monolithic_Solve_ComputeWetVolume")
        #WetVolumeBeforeRedistance = BiphasicFillingUtilities().ComputeWetVolume(self.model_part)
        #WetVolumeAfterConvecting=WetVolumeBeforeRedistance
        #Timer.Stop("DPG_Monolithic_Solve_ComputeWetVolume")
        # Compute the Flow Length
        # Checking the Conditions for the Redistance:CDL
        Timer.Start("DPG_Monolithic_Solve_DoRedistance")
        redistance_now=False
        if(self.internal_step_counter >= self.next_redistance):
            redistance_now=True
        else:
            min_open_node_distance=self.distance_utilities.CheckForRedistance(self.model_part)
            if (min_open_node_distance==0.0): #<(self.CFL*self.max_edge_size)):
                redistance_now=True

        if(redistance_now==True):
            #net_volume = self.model_part.ProcessInfo[NET_INPUT_MATERIAL]
            #ensure that inlet nodes are still wet
            for node in self.inlet_nodes:
                node.Free(DISTANCE)
                if( node.GetSolutionStepValue(DISTANCE) > 0.0):
                    node.SetSolutionStepValue(DISTANCE,0,  -0.1*self.max_edge_size)
            self.DoRedistance()
            self.next_redistance = self.internal_step_counter + self.redistance_frequency
        Timer.Stop("DPG_Monolithic_Solve_DoRedistance")

        # Now we compute the volume before the correction
        #Timer.Start("DPG_Monolithic_Solve_ComputeWetVolume")
        #WetVolumeBeforeCorrection = BiphasicFillingUtilities().ComputeWetVolume(self.model_part)
        #WetVolumeAfterRedistance=WetVolumeBeforeCorrection
        # Here The Net InleVolume is computed and saved into ProcessInfo[NET_INPUT_MATERIAL]
        #Timer.Stop("DPG_Monolithic_Solve_ComputeWetVolume")
        Timer.Start("DPG_Monolithic_Solve_CorrectVolume")
        self.CorrectVolume()
        #BiphasicFillingUtilities().ComputeNetInletVolume(self.model_part)
        #Timer.Start("DPG_Monolithic_Solve_CorrectVolume_OLD")
        #if(self.volume_correction_switch and self.step > self.vol_cr_step):
        #    net_volume = self.model_part.ProcessInfo[NET_INPUT_MATERIAL]
        #    self.volume_corrector.Execute(net_volume)
        #    #print("Computing Volume Correction:")
        #    #print("Average h: " + str(self.average_edge_size) )
        #    #print("Max h: " + str(self.max_edge_size) )
        #    #print("Net volume: "+str(net_volume))
        #    #BiphasicFillingUtilities().VolumeCorrection(self.model_part, net_volume, self.average_edge_size,self.negative_volume_correction)
        #else:
        #    self.model_part.ProcessInfo[WET_VOLUME] = BiphasicFillingUtilities().ComputeWetVolume(self.model_part) #fluid_model_part.ProcessInfo[WET_VOLUME]
        #Timer.Stop("DPG_Monolithic_Solve_CorrectVolume_OLD")
        Timer.Stop("DPG_Monolithic_Solve_CorrectVolume")

        Timer.Start("DPG_Monolithic_Solve_ApplyFluidProperties")
        self.ApplyFluidProperties()
        BiphasicFillingUtilities().ViscosityBasedSolidification(self.model_part,100.0)
        #self.IncreaseCSmagToSOlidify(50.0)
        Timer.Stop("DPG_Monolithic_Solve_ApplyFluidProperties")

        Timer.Start("DPG_Monolithic_Solve_Solve")
        #if(step > 3):
        #    (self.solver).Predict()
            #self.velocity_prediction.PredictVelocity()
        #ActivationUtilities().ActivateElementsAndConditions( self.model_part, DISTANCE, self.max_distance, True)
        print("")
        print("* Solving Navier-Stokes")
        (self.solver).Solve()

        self.internal_step_counter += 1
        Timer.Stop("DPG_Monolithic_Solve_Solve")

        # Now we compute the volume before the correction
        Timer.Start("DPG_Monolithic_Solve_ComputeMassLoss")
        WetVolumeAfterSolution = BiphasicFillingUtilities().ComputeWetVolume(self.model_part)
        #DischargeLossInWalls,DischargeLossInWetWalls=self.ComputeMassLossInWalls(self.model_part) #pOSITIVE OUTWARD, MASS LOSS. THE INPUT IS NEGATIVE
        #Dt=self.model_part.ProcessInfo[DELTA_TIME]
        #log_file = open("volume_loss.log", 'a')
        #volumen_inyectado=self.model_part.ProcessInfo[NET_INPUT_MATERIAL]-self.OldNetInletVolume
        #self.AcumulatedMassLossInWetWalls+=(DischargeLossInWetWalls*Dt)#+volumen_inyectado #wE HAVE TO CORRECT THE INPUT MATERIAL
        #VolTeor=WetVolumeBeforeConvecting+volumen_inyectado

        #self.AcumulatedConvectedMass+=( WetVolumeAfterConvecting -WetVolumeBeforeConvecting)
        #self.AcumulatedLossInRedistance+=(WetVolumeAfterRedistance-WetVolumeBeforeRedistance)

        #print("B. Convec: ", WetVolumeBeforeConvecting, " A.Conv: " , WetVolumeBeforeRedistance," Teor: ",VolTeor,"  A. Redist: ", WetVolumeBeforeCorrection)
        #log_file.write("t: "+str(self.model_part.ProcessInfo[TIME])+" Injected: " + str(-1.0*self.model_part.ProcessInfo[NET_INPUT_MATERIAL]) +" B.Convec: "+ str(WetVolumeBeforeConvecting) + " A.Conv: " + str(WetVolumeBeforeRedistance) +" Teor: "+str(VolTeor)+"  A.Redist: "+ str(WetVolumeBeforeCorrection)+" Loss.T.Walls: "+str(-1.0*self.AcumulatedMassLossInWetWalls)+"\n")
        #log_file.write("t: "+str(self.model_part.ProcessInfo[TIME])+" Injected: " + str(-1.0*self.model_part.ProcessInfo[NET_INPUT_MATERIAL])+ " Acum.Convected.Mass: "+str(self.AcumulatedConvectedMass)+ " A.ChangeRedistance: " + str(self.AcumulatedLossInRedistance)+" Loss.T.Walls: "+str(-1.0*self.AcumulatedMassLossInWetWalls)+" RealVolume: " +str(WetVolumeAfterSolution) +"\n")
        Timer.Stop("DPG_Monolithic_Solve_ComputeMassLoss")
        Timer.Stop("DPG_Monolithic_Solve")
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


    def ComputeMassLossInWalls(self, model_part):
        mass=0.0
        mass_on_fluid=0.0
        avg_vel=Vector(3)

        # Obtain mass through the outer surface for any condition
        # for condition in model_part.Conditions:
            # avg_vel[0]=0.0
            # avg_vel[1]=0.0
            # avg_vel[2]=0.0
            # area_normal=condition.GetValue(NORMAL)
            # area_norm=math.sqrt(area_normal*area_normal)
            # is_condition_wet=True
            # for node in condition.GetNodes():
                # tmp=node.GetSolutionStepValue(VELOCITY)
                # avg_vel=avg_vel+tmp
                # if(node.GetSolutionStepValue(DISTANCE)>=0):
                    # is_condition_wet=False
            # avg_vel=avg_vel*(1.0/3.0)
            # q=avg_vel*area_normal

            # mass+=q
            # if(is_condition_wet==True):
                # mass_on_fluid+=q
        return mass,mass_on_fluid
