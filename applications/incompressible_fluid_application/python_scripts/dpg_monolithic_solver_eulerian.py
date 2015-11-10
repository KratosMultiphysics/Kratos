from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.IncompressibleFluidApplication import *
from KratosMultiphysics.ThermoMechanicalApplication import *
from KratosMultiphysics.MeshingApplication import *
from KratosMultiphysics.FluidDynamicsApplication import *
from KratosMultiphysics.ExternalSolversApplication import *
from KratosMultiphysics.Click2CastApplication import *
# Check that KratosMultiphysics was imported in the main script
CheckForPreviousImport()

#import levelset_solver
import linear_solver_factory

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
    #model_part.AddNodalSolutionStepVariable(NODAL_MASS_BALANCE) #
    # variables needed for the distance solver
    #levelset_solver.AddVariables(model_part, distance_settings)
    model_part.AddNodalSolutionStepVariable(DISTANCE)

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
    #levelset_solver.AddDofs(model_part, distance_settings)
    print("dofs for the monolithic solver added correctly")


class MonolithicSolver:

    def __init__(self, model_part, domain_size,fluid_linear_solver_settings,redistance_settings):
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
        self.gmres_size = fluid_linear_solver_settings.gmres_size
        self.iterations = fluid_linear_solver_settings.linear_solver_iterations #400 # Ojo, antes 200
        self.tol = fluid_linear_solver_settings.linear_solver_tolerance #1e-5 #Before 1e-5
        self.verbosity = fluid_linear_solver_settings.verbosity
        self.linear_solver = ScalingSolver(AMGCLSolver(
            AMGCLSmoother.ILU0,
            AMGCLIterativeSolverType.BICGSTAB_WITH_GMRES_FALLBACK, #AMGCLIterativeSolverType.GMRES, #BICGSTAB_WITH_GMRES_FALLBACK,
            fluid_linear_solver_settings.linear_solver_tolerance,
            fluid_linear_solver_settings.linear_solver_iterations,
            fluid_linear_solver_settings.verbosity,
            fluid_linear_solver_settings.gmres_size),True)
        print("#####################################")
        print("#### LINEAR SOLVER               ####")
        print("#####################################")
        print(self.linear_solver)
        print("#####################################")
#### PRUEBA a ver 

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

# print "Construction monolithic solver finished"

        # Creat Lavel_set solver
        # construct the model part
        #if(domain_size == 2):
            #raise "error, still not implemented in 2D"
            #conv_elem = "SUPGConv2D"
            #conv_cond = "Condition2D"
        #else:
            #conv_elem = "SUPGConv3D" #"SUPGConvLevelSet"#"SUPGConv3D"
            #conv_cond = "Condition3D"
        #self.level_set_model_part = ModelPart("level_set_model_part")
        #self.conv_generator = ConnectivityPreserveModeler()
        #(self.conv_generator).GenerateModelPart(self.model_part,self.level_set_model_part, conv_elem, conv_cond)
        ##(ParallelFillCommunicator(self.level_set_model_part)).Execute();

        ## constructing the convection solver for the distance
        #self.level_set_solver = levelset_solver.Solver(self.level_set_model_part,domain_size,distance_settings)
        #self.level_set_solver.max_iter = 8
        
        max_cfl = 3.0;
        linear_solver = ScalingSolver(AMGCLSolver(
            AMGCLSmoother.ILU0,
            AMGCLIterativeSolverType.GMRES, #AMGCLIterativeSolverType.GMRES, #BICGSTAB_WITH_GMRES_FALLBACK,
            1e-9,
            50,
            fluid_linear_solver_settings.verbosity,
            50),True)
        #linear_solver = BICGSTABSolver(1e-9,200,DiagonalPreconditioner())
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
        # self.max_distance = self.max_edge_size * self.redistance_settings.max_distance_factor #Ojo antes 5.0



        ##########################################
        ##### Compute Max Edge Size       ########
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


            self.max_edge_size = self.distance_calculator.FindMaximumEdgeSize(self.level_set_model_part)
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

        # element size
        self.maxmin = []
        ParticleLevelSetUtils3D().FindMaxMinEdgeSize(
            self.model_part, self.maxmin)
        # Variables needed for computing the Efficiency of the Injection
        self.OldNetInletVolume=0.0
        self.OldWetVolume=0.0

        self.distance_utilities=DistanceUtilities()
    #

    def ApplyFluidProperties(self):
        # apply density
        mu1 = 1.0 * self.mu / self.rho1
        # mu1 = self.mu
        # mu2 = 0.01*self.mu/self.rho2
        mu2 = mu1
        BiphasicFillingUtilities().ApplyFluidProperties(self.model_part, mu1, self.rho1, mu2, self.rho2)

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
        print(">>>>>>>>>>>>>>> OSS_SWITCH = ", self.oss_switch)
        print("----------------------------- dynamic tau fluid -------------------------------------")
        self.model_part.ProcessInfo.SetValue(
            DYNAMIC_TAU, self.dynamic_tau_fluid)
        self.model_part.ProcessInfo.SetValue(OSS_SWITCH, self.oss_switch)

        # LEvel_set solver initialization
        #self.level_set_solver.dynamic_tau = self.dynamic_tau_levelset

        #Calculate Distances

        #self.redistance_utils.CalculateDistances(self.model_part,DISTANCE, NODAL_AREA,self.max_levels,self.max_distance)
        self.DoRedistance()

        #self.level_set_solver.linear_solver = AMGCLSolver(
            #AMGCLSmoother.ILU0,
            #AMGCLIterativeSolverType.GMRES,
            #1e-6, #self.tol,
            #200,
            #self.verbosity,
            #self.gmres_size)
        #self.level_set_solver.Initialize()

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
            #if(self.swap_normals==True):
            #    print(".............................")
            #    print("...Performing Normal Swapping")
            #    print(".............................")
            #    self.normal_util.SwapNormals(self.model_part)
            # Now we compute Normals
            self.normal_util.CalculateOnSimplex(self.model_part,self.domain_size, IS_STRUCTURE, 0, 35.0)  # ,0.0),180) #35.0)  # ,0.0,35.0

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
        if(self.redistance_type=="Old"):
            print("performing Old Redistance")
            self.distance_calculator.CalculateDistances(self.model_part,DISTANCE,NODAL_AREA, self.max_levels, self.max_distance)
        else:
            print("performing New Redistance")
            self.distance_calculator.Execute()

    def ConvectDistance(self):
        #self.level_set_model_part.ProcessInfo = self.model_part.ProcessInfo
        #(self.level_set_model_part.ProcessInfo).SetValue(CONVECTION_DIFFUSION_SETTINGS, distance_settings)
        #(self.level_set_model_part.ProcessInfo).SetValue(DYNAMIC_TAU, self.dynamic_tau_levelset)  # self.dynamic_tau
        #(self.level_set_solver).Solve()
        #BiphasicFillingUtilities().DistanceFarRegionCorrection(self.model_part,self.max_distance)
        
        self.levelset_convection_process.Execute()
        BiphasicFillingUtilities().DistanceFarRegionCorrection(self.model_part,self.max_distance)

    def Solve(self, step):
        # at the beginning of the calculations do a div clearance step
        if(self.divergence_clearance_performed == False):
            for node in self.model_part.Nodes:
                node.SetSolutionStepValue(DISTANCE,1,node.GetSolutionStepValue(DISTANCE))
            self.divergence_clearance_performed = True

        # Now we store the Old Values
        self.OldNetInletVolume=self.model_part.ProcessInfo[NET_INPUT_MATERIAL]
        self.OldWetVolume=self.model_part.ProcessInfo[WET_VOLUME]

        WetVolumeBeforeConvecting=BiphasicFillingUtilities().ComputeWetVolume(self.model_part)
        if(step > 3):
            (self.solver).Predict()

        Timer.Start("ConvectDistance")
        # convect distance function
        
        self.ConvectDistance()
        # recompute distance function as needed
        Timer.Start("DoRedistance")
        # Time to compute the Element Distance Gradient
        #self.distance_utilities.ComputeElementalGradient(self.model_part)
        WetVolumeBeforeRedistance = BiphasicFillingUtilities().ComputeWetVolume(self.model_part)
        # Checking the Conditions for the Redistance:CDL
        redistance_now=False
        if(self.internal_step_counter >= self.next_redistance):
            redistance_now=True
            print("")
            print("")
            print("Forced Redistance due to the number os time steps without redistance")
            print("")
            print("")
        else:
            min_open_node_distance=self.distance_utilities.CheckForRedistance(self.model_part)
            print(" min_open_node_distance %f" %min_open_node_distance)
            if (min_open_node_distance==0.0): #<(self.CFL*self.max_edge_size)):
                redistance_now=True
                print("Forced Redistance CFL=%f" % self.CFL, " Max_Edge= %f " %self.max_edge_size, " Min_open_node_distance= %f" %min_open_node_distance)

        if(redistance_now==True):
            #net_volume = self.model_part.ProcessInfo[NET_INPUT_MATERIAL]
            #BiphasicFillingUtilities().VolumeCorrection(self.model_part, net_volume, self.max_edge_size)
            
            #ensure that inlet nodes are still wet
            for node in self.inlet_nodes:
                node.Free(DISTANCE)
                if( node.GetSolutionStepValue(DISTANCE) > 0.0):
                    node.SetSolutionStepValue(DISTANCE,0,  -0.01*self.max_edge_size)
                    #pass
            self.DoRedistance()
            self.next_redistance = self.internal_step_counter + self.redistance_frequency
        
        Timer.Stop("DoRedistance")

        # Now we compute the volume before the correction
        WetVolumeBeforeCorrection = BiphasicFillingUtilities().ComputeWetVolume(self.model_part)
        
        # Here The Net InleVolume is computed and saved into ProcessInfo[NET_INPUT_MATERIAL]
        BiphasicFillingUtilities().ComputeNetInletVolume(self.model_part)
        log_file = open("volume_loss.log", 'a')
        volumen_inyectado=self.model_part.ProcessInfo[NET_INPUT_MATERIAL]-self.OldNetInletVolume
        VolTeor=WetVolumeBeforeConvecting+volumen_inyectado
        print("B. Convec: ", WetVolumeBeforeConvecting, " A. Conv: " , WetVolumeBeforeRedistance," Teor: ",VolTeor,"  A. Redist: ", WetVolumeBeforeCorrection)
        log_file.write("t: "+str(self.model_part.ProcessInfo[TIME]) +" B. Convec: "+ str(WetVolumeBeforeConvecting) + " A. Conv: " + str(WetVolumeBeforeRedistance) +" Teor: "+str(VolTeor)+"  A. Redist: "+ str(WetVolumeBeforeCorrection)+"\n")
        # Here The Net InleVolume is computed and saved into ProcessInfo[NET_INPUT_MATERIAL]
        # Now we compute the efficiency
        Efficiency=(WetVolumeBeforeCorrection-self.OldWetVolume)/(self.model_part.ProcessInfo[NET_INPUT_MATERIAL]-self.OldNetInletVolume)
        Efficiency=max(0,min(Efficiency,1.0))
        self.model_part.ProcessInfo[INJECTION_EFFICIENCY]=Efficiency


        if(self.volume_correction_switch and step > self.vol_cr_step):
            net_volume = self.model_part.ProcessInfo[NET_INPUT_MATERIAL]
            BiphasicFillingUtilities().VolumeCorrection(self.model_part, net_volume, self.max_edge_size,self.negative_volume_correction)
            print("Performed Volume Correction")
        else: # Ojo que lo acabo de
            self.model_part.ProcessInfo[WET_VOLUME] = BiphasicFillingUtilities().ComputeWetVolume(self.model_part) #fluid_model_part.ProcessInfo[WET_VOLUME]

        
        Timer.Start("ApplyFluidProperties")
        self.ApplyFluidProperties()
        BiphasicFillingUtilities().ViscosityBasedSolidification(self.model_part,100.0)
        #self.IncreaseCSmagToSOlidify(50.0)
        Timer.Stop("ApplyFluidProperties")

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
