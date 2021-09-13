from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

#Import basic python Libraries
import numpy as np
import math

# Timer
import time

# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.kratos_utilities as KratosUtilities

# Import applications
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD
have_conv_diff = KratosUtilities.CheckIfApplicationsAvailable("ConvectionDiffusionApplication")
if have_conv_diff:
    import KratosMultiphysics.ConvectionDiffusionApplication as KratosConvDiff

# Import base class file
from KratosMultiphysics.FluidDynamicsApplication.fluid_solver import FluidSolver
from KratosMultiphysics.FluidDynamicsApplication.read_distance_from_file import DistanceImportUtility

import KratosMultiphysics.python_linear_solver_factory as linear_solver_factory

def CreateSolver(model, custom_settings):
    return NavierStokesTwoFluidsSolver(model, custom_settings)

class NavierStokesTwoFluidsSolver(FluidSolver):

    @classmethod
    def GetDefaultSettings(cls):
        ##settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
        {
            "solver_type": "two_fluids_solver_from_defaults",
            "model_part_name": "",
            "domain_size": -1,
            "model_import_settings": {
                "input_type": "mdpa",
                "input_filename": "unknown_name",
                "reorder": false
            },
            "material_import_settings": {
                "materials_filename": ""
            },
            "distance_reading_settings"    : {
                "import_mode"         : "from_mdpa",
                "distance_file_name"  : "no_distance_file"
            },
            "maximum_iterations": 7,
            "echo_level": 0,
            "time_order": 2,
            "compute_reactions": false,
            "reform_dofs_at_each_step": false,
            "relative_velocity_tolerance": 1e-3,
            "absolute_velocity_tolerance": 1e-5,
            "relative_pressure_tolerance": 1e-3,
            "absolute_pressure_tolerance": 1e-5,
            "linear_solver_settings"       : {
                "solver_type"         : "amgcl"
            },
            "volume_model_part_name" : "volume_model_part",
            "skin_parts": [""],
            "assign_neighbour_elements_to_conditions": true,
            "no_skin_parts":[""],
            "time_stepping"                : {
                "automatic_time_step" : true,
                "CFL_number"          : 1,
                "minimum_delta_time"  : 1e-2,
                "maximum_delta_time"  : 1.0,
                "time_step"           : 0.0
            },
            "periodic": "periodic",
            "move_mesh_flag": false,
            "formulation": {
                "dynamic_tau": 1.0
            },
            "bfecc_convection" : false,
            "bfecc_number_substeps" : 10
        }""")

        default_settings.AddMissingParameters(super(NavierStokesTwoFluidsSolver, cls).GetDefaultSettings())
        return default_settings

    def __init__(self, model, custom_settings):
        self._validate_settings_in_baseclass=True # To be removed eventually
        super(NavierStokesTwoFluidsSolver,self).__init__(model,custom_settings)

        self.element_name = "TwoFluidNavierStokes"
        self.condition_name = "NavierStokesWallCondition"
        self.element_has_nodal_properties = True

        self.min_buffer_size = 3

        self._bfecc_convection = self.settings["bfecc_convection"].GetBool()

        ## Set the distance reading filename
        # TODO: remove the manual "distance_file_name" set as soon as the problem type one has been tested.
        if (self.settings["distance_reading_settings"]["import_mode"].GetString() == "from_GiD_file"):
            self.settings["distance_reading_settings"]["distance_file_name"].SetString(self.settings["model_import_settings"]["input_filename"].GetString()+".post.res")

        KratosMultiphysics.Logger.PrintInfo("NavierStokesTwoFluidsSolver", "Construction of NavierStokesTwoFluidsSolver finished.")

    def AddVariables(self):
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DENSITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DYNAMIC_VISCOSITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PRESSURE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ACCELERATION)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.MESH_VELOCITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.IS_STRUCTURE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.BODY_FORCE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_H)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_AREA)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION_WATER_PRESSURE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.EXTERNAL_PRESSURE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.FLAG_VARIABLE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE)              # Distance function nodal values
        self.main_model_part.AddNodalSolutionStepVariable(KratosCFD.DISTANCE_AUX)                   # Auxiliary distance function nodal values
        self.main_model_part.AddNodalSolutionStepVariable(KratosCFD.DISTANCE_AUX2)                  # Auxiliary distance function nodal values
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE_GRADIENT)     # Distance gradient nodal values
        self.main_model_part.AddNodalSolutionStepVariable(KratosCFD.DISTANCE_GRADIENT_AUX)          # Auxiliary Distance gradient nodal values
        self.main_model_part.AddNodalSolutionStepVariable(KratosCFD.CONVECTIVE_VELOCITY)            # Store conctive velocity for level-set process
        self.main_model_part.AddNodalSolutionStepVariable(KratosCFD.CURVATURE)                      # Store curvature as a nodal variable
        self.main_model_part.AddNodalSolutionStepVariable(KratosCFD.AREA_VARIABLE_AUX)              # Auxiliary area_variable for parallel distance calculator    
        self.main_model_part.AddNodalSolutionStepVariable(KratosCFD.NORMAL_VECTOR)                  # Auxiliary normal vector at interface
        self.main_model_part.AddNodalSolutionStepVariable(KratosCFD.TANGENT_VECTOR)                 # Auxiliary tangent vector at contact line
        self.main_model_part.AddNodalSolutionStepVariable(KratosCFD.CONTACT_VECTOR)                 # Auxiliary contact vector 
        self.main_model_part.AddNodalSolutionStepVariable(KratosCFD.CONTACT_ANGLE)                  # Contact angle (may not be needed at nodes)  
        self.main_model_part.AddNodalSolutionStepVariable(KratosCFD.CONTACT_VECTOR_MICRO)           # Auxiliary contact vector at micro-scale
        self.main_model_part.AddNodalSolutionStepVariable(KratosCFD.CONTACT_ANGLE_MICRO)            # Contact angle (micro-scale)
        self.main_model_part.AddNodalSolutionStepVariable(KratosCFD.CONTACT_VELOCITY)               # Contact line tangential velocity (normal to the contact-line) 
        #self.main_model_part.AddNodalSolutionStepVariable(KratosCFD.VELOCITY_STAR)                  # Last known velocity
        #self.main_model_part.AddNodalSolutionStepVariable(KratosCFD.PRESSURE_STAR)                  # Last known pressure
        self.main_model_part.AddNodalSolutionStepVariable(KratosCFD.PRESSURE_GRADIENT_AUX)          # Pressure gradient on positive and negative sides

        KratosMultiphysics.Logger.PrintInfo("NavierStokesTwoFluidsSolver", "Fluid solver variables added correctly.")

    def PrepareModelPart(self):
        # Initialize the level-set function
        if not self.main_model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED]:
            ## Setting the nodal distance
            self._set_distance_function()

        # Call the base solver PrepareModelPart()
        super(NavierStokesTwoFluidsSolver, self).PrepareModelPart()

    def Initialize(self):
        self.computing_model_part = self.GetComputingModelPart()

        # Geting the slip condition model part
        if self.model.HasModelPart("FluidModelPart.Slip3D"):
            self.slip_model_part = self.model.GetModelPart("FluidModelPart.Slip3D")
        #     for slip_condition in self.slip_model_part.Conditions:
        #         for node in slip_condition.GetNodes():
        #             print(node.X)

        # Geting the injection condition model part
        if self.model.HasModelPart("FluidModelPart.AutomaticInlet3D_Injection"):
            self.injection_model_part = self.model.GetModelPart("FluidModelPart.AutomaticInlet3D_Injection")
            for injection_condition in self.injection_model_part.Conditions:
                injection_condition.Set(KratosMultiphysics.BOUNDARY, True)
                for node in injection_condition.GetNodes():
                    node.Set(KratosMultiphysics.BOUNDARY, True)
                    if node.GetSolutionStepValue(KratosMultiphysics.DISTANCE) > -1.0e-6:
                        node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, -1.0e-6)

        ## Construct the linear solver
        self.linear_solver = linear_solver_factory.ConstructSolver(self.settings["linear_solver_settings"])

        KratosMultiphysics.NormalCalculationUtils().CalculateOnSimplex(self.computing_model_part, self.computing_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE])

        self.neighbour_search = KratosMultiphysics.FindNodalNeighboursProcess(self.computing_model_part, 10, 10)
        (self.neighbour_search).Execute()

        self.accelerationLimitationUtility = KratosCFD.AccelerationLimitationUtilities( self.computing_model_part, 5.0 )

        # If needed, create the estimate time step utility
        if (self.settings["time_stepping"]["automatic_time_step"].GetBool()):
            self.EstimateDeltaTimeUtility = self._GetAutomaticTimeSteppingUtility()

        # Set the time discretization utility to compute the BDF coefficients
        time_order = self.settings["time_order"].GetInt()
        if time_order == 2:
            self.time_discretization = KratosMultiphysics.TimeDiscretization.BDF(time_order)
        else:
            raise Exception("Only \"time_order\" equal to 2 is supported. Provided \"time_order\": " + str(time_order))

        # Creating the solution strategy
        self.conv_criteria = KratosCFD.VelPrCriteria(self.settings["relative_velocity_tolerance"].GetDouble(),
                                                     self.settings["absolute_velocity_tolerance"].GetDouble(),
                                                     self.settings["relative_pressure_tolerance"].GetDouble(),
                                                     self.settings["absolute_pressure_tolerance"].GetDouble())

        (self.conv_criteria).SetEchoLevel(self.settings["echo_level"].GetInt())

        self.find_neighbouring_elements_process = self._set_find_neighbouring_elements_process()
        (self.find_neighbouring_elements_process).Execute()

        #self.find_neighbouring_nodes_process = self._set_find_neighbouring_nodes_process()
        #(self.find_neighbouring_nodes_process).Execute()

        self.level_set_convection_process = self._set_level_set_convection_process()

        #Set IS_STRUCTURE to define contact line: 0.0: not needed, and 1.0: moved to apply_slip_condition
        #for node in self.main_model_part.Nodes:
        #    node.SetValue(KratosMultiphysics.IS_STRUCTURE, 0.0)

        #self.mass_conservation_correction = self._set_mass_conservation_correction()
        #(self.mass_conservation_correction).Initialize();

        self.distance_gradient_process = self._set_distance_gradient_process()
        (self.distance_gradient_process).Execute()

        self.curvature_calculation_process = self._set_curvature_calculation_process()
        #(self.curvature_calculation_process).Execute()

        self.interface_curvature_calculation = self._set_interface_curvature_calculation()
        #(self.interface_curvature_calculation).Execute()

        KratosMultiphysics.Logger.PrintInfo("NavierStokesTwoFluidsSolver", "Start re-distancing")

        self.hyperbolic_distance_reinitialization = self._set_hyperbolic_distance_reinitialization()
        #(self.hyperbolic_distance_reinitialization).Execute()

        self.parallel_distance_process = self._set_parallel_distance_process()
        #layers = int(2000/100000*self.main_model_part.NumberOfElements())
        #(self.parallel_distance_process).CalculateDistances(
        #            self.main_model_part,
        #            KratosMultiphysics.DISTANCE,
        #            KratosCFD.AREA_VARIABLE_AUX,
        #            layers,
        #            1.0e0,
        #            (self.parallel_distance_process).NOT_CALCULATE_EXACT_DISTANCES_TO_PLANE) #NOT added on feb 20, 2020

        #self.variational_distance_process = self._set_variational_distance_process()
        #(self.variational_distance_process).Execute()

        #self.lumped_eikonal_distance_calculation = self._set_lumped_eikonal_distance_calculation()
        #(self.lumped_eikonal_distance_calculation).Execute()

        #self.variational_non_eikonal_distance = self._set_variational_non_eikonal_distance()
        #(self.distance_gradient_process).Execute()
        #(self.curvature_calculation_process).Execute()
        #(self.variational_non_eikonal_distance).Execute()
        #for node in self.main_model_part.Nodes:
        #    smooth_distance = node.GetSolutionStepValue(KratosCFD.DISTANCE_AUX2)
        #    node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, smooth_distance)

        #KratosMultiphysics.Logger.PrintInfo("NavierStokesTwoFluidsSolver", "Re-distancing is finished")

        self.surface_smoothing_process = self._set_surface_smoothing_process()

        print("Contact Angle Evaluator Const 1")
        self.contact_angle_evaluator = self._set_contact_angle_evaluator()
        print("Contact Angle Evaluator Const 2")

        #print("Smoothing")
        #print(time.time())
        #(self.surface_smoothing_process).Execute()
        #print(time.time())
        #for node in self.main_model_part.Nodes:
        #    smooth_distance = node.GetSolutionStepValue(KratosCFD.DISTANCE_AUX)
        #    node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, smooth_distance)

        #it_number=self.linear_solver.GetIterationsNumber()
        #KratosMultiphysics.Logger.PrintInfo("linear solver number of iterations, smoothing", it_number)

        self.interface_pressure_gradient_calculation = self._set_interface_pressure_gradient_calculation()
        #(self.interface_pressure_gradient_calculation).Execute()

        time_scheme = KratosMultiphysics.ResidualBasedIncrementalUpdateStaticSchemeSlip(self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE],   # Domain size (2,3)
                                                                                        self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]+1) # DOFs (3,4)

        builder_and_solver = KratosMultiphysics.ResidualBasedBlockBuilderAndSolver(self.linear_solver)

        self.solver = KratosMultiphysics.ResidualBasedNewtonRaphsonStrategy(self.computing_model_part,
                                                                            time_scheme,
                                                                            self.linear_solver,
                                                                            self.conv_criteria,
                                                                            builder_and_solver,
                                                                            self.settings["maximum_iterations"].GetInt(),
                                                                            self.settings["compute_reactions"].GetBool(),
                                                                            self.settings["reform_dofs_at_each_step"].GetBool(),
                                                                            self.settings["move_mesh_flag"].GetBool())

        (self.solver).SetEchoLevel(self.settings["echo_level"].GetInt())

        (self.solver).Initialize() # Initialize the solver. Otherwise the constitutive law is not initializated.
        (self.solver).Check()

        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DYNAMIC_TAU, self.settings["formulation"]["dynamic_tau"].GetDouble())

        with open("ZeroDistance_Structured.log", "w") as SdistLogFile:
            SdistLogFile.write( "time_step" + "\t" + "XZeroMin" + "\t" + "XZeroMax" + "\t" + "ZZero" + "\n" )

        with open("ZeroDistance_Unstructured.log", "w") as USdistLogFile:
            USdistLogFile.write( "time_step" + "\t" + "MeanRadius" + "\n" )

        with open("solver_iteration.log", "w") as iterLogFile:
            iterLogFile.write( "time_step" + "\t" + "iter_number" + "\n" )

        with open("MaxVelocity.log", "w") as velLogFile:
            velLogFile.write( "time_step" + "\t" + "VMax" + "\n" )

        with open("ContactAngle.log", "w") as CangleLogFile:
            CangleLogFile.write( "element_id" + "\t" + "contact_angle" + "\n" )

        with open("MeanContactAngle.log", "w") as CangleLogFile:
            CangleLogFile.write( "time_step" + "\t" + "mean_contact_angle" + "\n" )

        with open("MeanContactAngleMicro.log", "w") as CangleLogFile:
            CangleLogFile.write( "time_step" + "\t" + "mean_contact_angle_micro" + "\n" )

        with open("ContactVelocity.log", "w") as CvelLogFile:
            CvelLogFile.write( "element_id" + "\t" + "contact_velocity" + "\n" )

        with open("MeanContactVelocity.log", "w") as CvelLogFile:
            CvelLogFile.write( "time_step" + "\t" + "mean_velocity" + "\n" )

        KratosMultiphysics.Logger.PrintInfo("NavierStokesTwoFluidsSolver", "Solver initialization finished.")

    def InitializeSolutionStep(self):
        KratosMultiphysics.Logger.PrintInfo("Navier Stokes Two Fluid Solver", "to begin InitializeSolutionStep")

        if self._TimeBufferIsInitialized():
            # (self.accelerationLimitationUtility).Execute()

            # Recompute the BDF2 coefficients
            (self.time_discretization).ComputeAndSaveBDFCoefficients(self.GetComputingModelPart().ProcessInfo)

            TimeStep = self.main_model_part.ProcessInfo[KratosMultiphysics.STEP]
            DT = self.main_model_part.ProcessInfo[KratosMultiphysics.DELTA_TIME]

            gravity = 9.81#0.0#TimeStep*DT/0.01*9.81
            if gravity > 9.81:
                gravity = 9.81

            tilting_angle = 0.0

            # if gravity == 9.81:
            #     tilting_angle = (TimeStep*DT - 0.01)/0.09*(90.0/180.0)*math.pi
            #     if tilting_angle < 0.0:
            #         tilting_angle = 0.0
            #     elif tilting_angle > (30.0/180.0)*math.pi:
            #         tilting_angle = (30.0/180.0)*math.pi

            sinAlpha = math.sin(tilting_angle)
            cosAlpha = math.cos(tilting_angle)

            for node in self.main_model_part.Nodes:
                node.SetSolutionStepValue(KratosMultiphysics.BODY_FORCE_X, gravity*sinAlpha)
                node.SetSolutionStepValue(KratosMultiphysics.BODY_FORCE_Y, 0.0)
                node.SetSolutionStepValue(KratosMultiphysics.BODY_FORCE_Z, -gravity*cosAlpha)

            # Recompute the distance field according to the new level-set position
            if (TimeStep % 1 == 0):
                print("Redistancing")
                print(time.time())

                #(self.variational_distance_process).Execute()

                #(self.hyperbolic_distance_reinitialization).Execute()

                layers = 200#int(4000/100000*self.main_model_part.NumberOfElements())
                (self.parallel_distance_process).CalculateInterfacePreservingDistances( #CalculateDistances(
                    self.main_model_part,
                    KratosMultiphysics.DISTANCE,
                    KratosCFD.AREA_VARIABLE_AUX,
                    layers,
                    1.0e0),
                    #(self.parallel_distance_process).CALCULATE_EXACT_DISTANCES_TO_PLANE)

                # print(time.time())

            #########################################
            ##
            ## Previous position of smoothing process
            ##
            #########################################

            print("Level-set")
            print(time.time())

            #if self.model.HasModelPart("FluidModelPart.AutomaticInlet3D_Injection"):
            #    for injection_condition in self.injection_model_part.Conditions:
            #        for node in injection_condition.GetNodes():
            #            if node.GetSolutionStepValue(KratosMultiphysics.DISTANCE) > -1.0e-6:
            #                node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, -1.0e-6)

            # Perform the level-set convection according to the previous step velocity
            if self._bfecc_convection:
                KratosMultiphysics.Logger.PrintInfo("LevelSetSolver", "BFECCconvect will be called")
                (self.level_set_convection_process).CopyScalarVarToPreviousTimeStep(
                    self.main_model_part,
                    KratosMultiphysics.DISTANCE)
                (self.level_set_convection_process).BFECCconvect(
                    self.main_model_part,
                    KratosMultiphysics.DISTANCE,
                    KratosMultiphysics.VELOCITY,
                    self.settings["bfecc_number_substeps"].GetInt())
            else:
                for node in self.main_model_part.Nodes:
                    velocityOld = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY, 1)
                    velocity = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY)
                    node.SetSolutionStepValue(KratosCFD.CONVECTIVE_VELOCITY, 1,
                        velocity )
                        #0.5 * (velocity + velocityOld) )
                        #velocityOld )
                    node.SetSolutionStepValue(KratosCFD.CONVECTIVE_VELOCITY,
                        velocity )
                        #0.5 * (velocity + velocityOld) )

                (self.level_set_convection_process).Execute()

            print(time.time())

            ####################################
            ##
            ## Original position of smoothing process
            ##
            ####################################
            # Smoothing the surface to filter oscillatory surface
            (self.distance_gradient_process).Execute() # Always check if calculated above
            print("Smoothing Started")
            print(time.time())
            (self.surface_smoothing_process).Execute()
            print(time.time())
            print("Smoothing Finished")

            TimeStep = self.main_model_part.ProcessInfo[KratosMultiphysics.STEP]
            if (TimeStep % 1 == 0):
                for node in self.main_model_part.Nodes:
                    smooth_distance = node.GetSolutionStepValue(KratosCFD.DISTANCE_AUX)
                    node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, smooth_distance)
                    #print(node.GetSolutionStepValue(KratosMultiphysics.BODY_FORCE))

            # Compute the DISTANCE_GRADIENT on nodes
            (self.distance_gradient_process).Execute()

            # Compute CURVATURE on nodes
            (self.curvature_calculation_process).Execute()

            # Evaluating the average nodal contact angle
            print("Contact Angle Evaluator Exe 1")
            (self.contact_angle_evaluator).Execute()
            print("Contact Angle Evaluator Exe 2")

            # Calculate nodal pressure gradient excluding the cut elements needed for stabilization of Kee_inv
            (self.interface_pressure_gradient_calculation).Execute()

            # Update the DENSITY and DYNAMIC_VISCOSITY values according to the new level-set
            self._SetNodalProperties()

            # Initialize the solver current step
            (self.solver).InitializeSolutionStep()

        KratosMultiphysics.Logger.PrintInfo("Navier Stokes Two Fluid Solver", "ended InitializeSolutionStep")

    def SolveSolutionStep(self):
        KratosMultiphysics.Logger.PrintInfo("Navier Stokes Two Fluid Solver", "to begin SolutionStep")

        # Do the necessary correction
        for node in self.main_model_part.Nodes:
            #node.Free(KratosMultiphysics.DISTANCE)
            dist = node.GetSolutionStepValue(KratosMultiphysics.DISTANCE)
            if (abs(dist) < 1.0e-12):
                print("Solver: do the correction")
                auxdist = 1.0e-5*node.GetValue(KratosMultiphysics.NODAL_H)
                print(auxdist)
                if (dist > 0.0):
                    node.SetSolutionStepValue(KratosMultiphysics.DISTANCE,auxdist)
                else:
                    node.SetSolutionStepValue(KratosMultiphysics.DISTANCE,-auxdist)

        # Call the base solver SolveSolutionStep()
        print("Navier-Stokes")
        print(time.time())
        super(NavierStokesTwoFluidsSolver, self).SolveSolutionStep()
        print(time.time())

        it_number=self.linear_solver.GetIterationsNumber()
        KratosMultiphysics.Logger.PrintInfo("linear solver number of iterations, NS", it_number)
        KratosMultiphysics.Logger.PrintInfo("Navier Stokes Two Fluid Solver", "solved SolutionStep")

        ####################################
        ##
        ## Most recent position of smoothing process
        ##
        ####################################

        it_number=self.linear_solver.GetIterationsNumber()
        KratosMultiphysics.Logger.PrintInfo("linear solver number of iterations, smoothing", it_number)

        if not(self._bfecc_convection):
            if self.model.HasModelPart("FluidModelPart.AutomaticInlet3D_Injection"):
                for injection_condition in self.injection_model_part.Conditions:
                    for node in injection_condition.GetNodes():
                        if node.GetSolutionStepValue(KratosMultiphysics.DISTANCE) > -1.0e-6:
                            node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, -1.0e-6)

            for node in self.main_model_part.Nodes:
                velocityOld = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY, 1)
                velocity = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY)
                node.SetSolutionStepValue(KratosCFD.CONVECTIVE_VELOCITY, 1,
                    0.5 * (velocity + velocityOld) )
                node.SetSolutionStepValue(KratosCFD.CONVECTIVE_VELOCITY,
                    velocity )

            (self.level_set_convection_process).Execute()

        # # Smoothing the surface to filter oscillatory surface
        # (self.distance_gradient_process).Execute() # Always check if calculated above
        # print("Smoothing")
        # print(time.time())
        # (self.surface_smoothing_process).Execute()
        # print(time.time())

        # TimeStep = self.main_model_part.ProcessInfo[KratosMultiphysics.STEP]
        # if (TimeStep % 500 == 0):
        #     for node in self.main_model_part.Nodes:
        #         smooth_distance = node.GetSolutionStepValue(KratosCFD.DISTANCE_AUX)
        #         node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, smooth_distance)

        """ n_pos = 0
        n_neg = 0
        for elem in self.main_model_part.Elements:
            for node in elem.GetNodes():
                dist = node.GetSolutionStepValue(KratosMultiphysics.DISTANCE)
                if (dist > 0.0):
                    n_pos += 1
                else:
                    n_neg += 1

            if (n_pos > 0 and n_neg > 0):
                for node in elem.GetNodes():
                    smooth_distance = node.GetSolutionStepValue(KratosCFD.DISTANCE_AUX)
                    node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, smooth_distance) """

    def FinalizeSolutionStep(self):
        TimeStep = self.main_model_part.ProcessInfo[KratosMultiphysics.STEP]
        DT = self.main_model_part.ProcessInfo[KratosMultiphysics.DELTA_TIME]

        if (TimeStep % 10 == 0):

            ##############################
            # Contact Angle and Velocity #
            ##############################

            mean_Cangle = 0.0
            mean_Cangle_micro = 0.0
            mean_Cvel = 0.0
            num_C = 0
            with open("ContactAngle.log", "a") as CangleLogFile, open("ContactVelocity.log", "a") as CvelLogFile:
                CangleLogFile.write( "\n" + str(TimeStep*DT) + "\n" )
                CvelLogFile.write( "\n" + str(TimeStep*DT) + "\n" )
                for elem in self.main_model_part.Elements:
                    cangle = elem.GetValue(KratosCFD.CONTACT_ANGLE)
                    cangle_micro = elem.GetValue(KratosCFD.CONTACT_ANGLE_MICRO)
                    cvel = elem.GetValue(KratosCFD.CONTACT_VELOCITY)
                    if (cangle != 0.0):
                        CangleLogFile.write( str(elem.Id) + "\t" + str(cangle) + "\n" )
                        CvelLogFile.write( str(elem.Id) + "\t" + str(cvel) + "\n" )
                        mean_Cangle += cangle
                        mean_Cangle_micro += cangle_micro
                        mean_Cvel += cvel
                        num_C += 1

            if (num_C > 1):
                with open("MeanContactAngle.log", "a") as meanCangleLogFile, open("MeanContactAngleMicro.log", "a") as meanMicroCangleLogFile, open("MeanContactVelocity.log", "a") as meanCvelLogFile:
                    meanCangleLogFile.write( str(TimeStep*DT) + "\t" + str(mean_Cangle/num_C) + "\n" )
                    meanMicroCangleLogFile.write( str(TimeStep*DT) + "\t" + str(mean_Cangle_micro/num_C) + "\n" )
                    meanCvelLogFile.write( str(TimeStep*DT) + "\t" + str(mean_Cvel/num_C) + "\n" )

            ###############################
            # Zero Disance - Unstructured #
            ###############################

            if self.model.HasModelPart("FluidModelPart.Slip3D"):
                x_center = 0.5#4.0e-3#1.5e-3#
                y_center = x_center
                z_center = 0.0
                mean_radius = 0.0
                num_points = 0.0
                for slip_condition in self.slip_model_part.Conditions:
                    num_positive = 0
                    num_negative = 0
                    for node in slip_condition.GetNodes():
                        dist = node.GetSolutionStepValue(KratosMultiphysics.DISTANCE)
                        if (dist > 0.0):
                            num_positive += 1
                        else:
                            num_negative += 1

                    if ((num_positive == 2 and num_negative == 1) or (num_positive == 1 and num_negative == 2)):
                        x_positive = [-1.0, -1.0]
                        y_positive = [-1.0, -1.0]
                        z_positive = [-1.0, -1.0]
                        x_negative = [-1.0, -1.0]
                        y_negative = [-1.0, -1.0]
                        z_negative = [-1.0, -1.0]
                        abs_dist_positive = [-1.0, -1.0]
                        abs_dist_negative = [-1.0, -1.0]

                        positive_pos = 0
                        negative_pos = 0
                        for node in slip_condition.GetNodes():
                            dist = node.GetSolutionStepValue(KratosMultiphysics.DISTANCE)
                            if (dist > 0):
                                x_positive[positive_pos] = node.X
                                y_positive[positive_pos] = node.Y
                                z_positive[positive_pos] = node.Z
                                abs_dist_positive[positive_pos] = abs(dist)
                                positive_pos += 1
                            else:
                                x_negative[negative_pos] = node.X
                                y_negative[negative_pos] = node.Y
                                z_negative[negative_pos] = node.Z
                                abs_dist_negative[negative_pos] = abs(dist)
                                negative_pos += 1

                        x_zero = -1.0
                        y_zero = -1.0
                        z_zero = -1.0
                        if (positive_pos == 2):
                            x_zero = 0.5*((abs_dist_positive[0]*x_negative[0]+abs_dist_negative[0]*x_positive[0])/(abs_dist_positive[0]+abs_dist_negative[0])+(abs_dist_positive[1]*x_negative[0]+abs_dist_negative[0]*x_positive[1])/(abs_dist_positive[1]+abs_dist_negative[0]))
                            y_zero = 0.5*((abs_dist_positive[0]*y_negative[0]+abs_dist_negative[0]*y_positive[0])/(abs_dist_positive[0]+abs_dist_negative[0])+(abs_dist_positive[1]*y_negative[0]+abs_dist_negative[0]*y_positive[1])/(abs_dist_positive[1]+abs_dist_negative[0]))
                            z_zero = 0.5*((abs_dist_positive[0]*z_negative[0]+abs_dist_negative[0]*z_positive[0])/(abs_dist_positive[0]+abs_dist_negative[0])+(abs_dist_positive[1]*z_negative[0]+abs_dist_negative[0]*z_positive[1])/(abs_dist_positive[1]+abs_dist_negative[0]))
                        else:
                            x_zero = 0.5*((abs_dist_positive[0]*x_negative[0]+abs_dist_negative[0]*x_positive[0])/(abs_dist_positive[0]+abs_dist_negative[0])+(abs_dist_positive[0]*x_negative[1]+abs_dist_negative[1]*x_positive[0])/(abs_dist_positive[0]+abs_dist_negative[1]))
                            y_zero = 0.5*((abs_dist_positive[0]*y_negative[0]+abs_dist_negative[0]*y_positive[0])/(abs_dist_positive[0]+abs_dist_negative[0])+(abs_dist_positive[0]*y_negative[1]+abs_dist_negative[1]*y_positive[0])/(abs_dist_positive[0]+abs_dist_negative[1]))
                            z_zero = 0.5*((abs_dist_positive[0]*z_negative[0]+abs_dist_negative[0]*z_positive[0])/(abs_dist_positive[0]+abs_dist_negative[0])+(abs_dist_positive[0]*z_negative[1]+abs_dist_negative[1]*z_positive[0])/(abs_dist_positive[0]+abs_dist_negative[1]))

                        num_points += 1.0
                        mean_radius += math.sqrt((x_center-x_zero)**2+(y_center-y_zero)**2+(z_center-z_zero)**2)

                with open("ZeroDistance_Unstructured.log", "a") as USdistLogFile:
                    if (num_points > 0):
                        USdistLogFile.write( str(TimeStep*DT) + "\t" + str(mean_radius/num_points) + "\n" )
                    else:
                        USdistLogFile.write( str(TimeStep*DT) + "\t" + "0.0" + "\n" )

            #############################
            # Zero Disance - Structured #
            #############################

            X_c = 0.5#4.0e-3#1.5e-3#
            Z_max = 1.0#3.0e-3#2.0e-3#
            X_min = 0.0
            X_max = 2*X_c
            Y_c = X_c
            Z_c = X_c

            XPlusMin = X_c
            XMinusMin = X_min
            DistPlusMin = 1.0e5
            DistMinusMin = -1.0e5
            XZeroMin = (XPlusMin + XMinusMin) / 2.0

            XPlusMax = X_max
            XMinusMax = X_c
            DistPlusMax = 1.0e5
            DistMinusMax = -1.0e5
            XZeroMax = (XPlusMax + XMinusMax) / 2.0

            ZPlus = Z_max
            ZMinus = 0.0
            DistPlusZ = 1.0e5
            DistMinusZ = -1.0e5
            ZZero = (ZPlus + ZMinus) / 2.0

            for node in self.main_model_part.Nodes:
                NodeX = node.X
                NodeY = node.Y
                NodeZ = node.Z
                if (abs(NodeY - Y_c) < 1.0e-6 and abs(NodeZ - 0.0) < 1.0e-6):
                    if (NodeX <= X_c):
                        Dist = node.GetSolutionStepValue(KratosMultiphysics.DISTANCE)
                        if (Dist >= 0.0 and Dist <= DistPlusMin):
                            DistPlusMin = Dist
                            XPlusMin = NodeX
                        if (Dist <= 0.0 and Dist >= DistMinusMin):
                            DistMinusMin = Dist
                            XMinusMin = NodeX
                    if (NodeX >= X_c):
                        Dist = node.GetSolutionStepValue(KratosMultiphysics.DISTANCE)
                        if (Dist >= 0.0 and Dist <= DistPlusMax):
                            DistPlusMax = Dist
                            XPlusMax = NodeX
                        if (Dist <= 0.0 and Dist >= DistMinusMax):
                            DistMinusMax = Dist
                            XMinusMax = NodeX

                if (abs(NodeY - Y_c) < 1.0e-6 and abs(NodeX - X_c) < 1.0e-6):
                    if (NodeZ > Z_c): #condition only for oscillating droplet. Omit it for spreading droplet
                        Dist = node.GetSolutionStepValue(KratosMultiphysics.DISTANCE)
                        if (Dist >= 0.0 and Dist <= DistPlusZ):
                            DistPlusZ = Dist
                            ZPlus = NodeZ
                        if (Dist <= 0.0 and Dist >= DistMinusZ):
                            DistMinusZ = Dist
                            ZMinus = NodeZ

            if (abs(DistPlusMin - DistMinusMin) > 1.0e-15):
                XZeroMin = XMinusMin + (-DistMinusMin)/(DistPlusMin - DistMinusMin)*(XPlusMin - XMinusMin)
            else:
                XZeroMin = XMinusMin

            if (abs(DistPlusMax - DistMinusMax) > 1.0e-15):
                XZeroMax = XMinusMax + (-DistMinusMax)/(DistPlusMax - DistMinusMax)*(XPlusMax - XMinusMax)
            else:
                XZeroMax = XMinusMax

            if (abs(DistPlusZ - DistMinusZ) > 1.0e-15):
                ZZero = ZMinus + (-DistMinusZ)/(DistPlusZ - DistMinusZ)*(ZPlus - ZMinus)
            else:
                ZZero = ZMinus

            with open("ZeroDistance_Structured.log", "a") as distLogFile:
                distLogFile.write( str(TimeStep*DT) + "\t" + str(XZeroMin) + "\t" + str(XZeroMax) + "\t" + str(ZZero) + "\n" )

            ####################
            # Maximum Velocity #
            ####################

            VMax = 0.0

            for node in self.main_model_part.Nodes:
                VX = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_X)
                VY = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_Y)
                VZ = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_Z)
                normv = math.sqrt(VX**2 + VY**2 + VZ**2)

                if (normv > VMax):
                    VMax = normv

            with open("MaxVelocity.log", "a") as velLogFile:
                velLogFile.write( str(TimeStep*DT) + "\t" + str(VMax) + "\n" )

        if self._TimeBufferIsInitialized():
            (self.solver).FinalizeSolutionStep()
            if (TimeStep >= 3):
                KratosMultiphysics.Logger.PrintInfo("Navier Stokes Two Fluid Solver, TIMESTEP= ", TimeStep)
                (self.accelerationLimitationUtility).Execute()

    # TODO: Remove this method as soon as the subproperties are available
    def _SetPhysicalProperties(self):
        import os
        warn_msg  = '\nThe materials import mechanism used in the two fluids solver is DEPRECATED!\n'
        warn_msg += 'It will be removed to use the base fluid_solver.py one as soon as the subproperties are available.\n'
        KratosMultiphysics.Logger.PrintWarning('\n\x1b[1;31mDEPRECATION-WARNING\x1b[0m', warn_msg)

        # Check if the fluid properties are provided using a .json file
        materials_filename = self.settings["material_import_settings"]["materials_filename"].GetString()
        if (materials_filename != ""):
            with open(materials_filename,'r') as materials_file:
                materials = KratosMultiphysics.Parameters(materials_file.read())

            # Create and read an auxiliary materials file for each one of the fields
            for i_material in materials["properties"]:
                aux_materials = KratosMultiphysics.Parameters()
                aux_materials.AddEmptyArray("properties")
                aux_materials["properties"].Append(i_material)
                prop_id = i_material["properties_id"].GetInt()

                aux_materials_filename = materials_filename + "_" + str(prop_id) + ".json"
                with open(aux_materials_filename,'w') as aux_materials_file:
                    aux_materials_file.write(aux_materials.WriteJsonString())
                    aux_materials_file.close()

                aux_material_settings = KratosMultiphysics.Parameters("""{"Parameters": {"materials_filename": ""}} """)
                aux_material_settings["Parameters"]["materials_filename"].SetString(aux_materials_filename)
                KratosMultiphysics.ReadMaterialsUtility(aux_material_settings, self.model)

                os.remove(aux_materials_filename)

            materials_imported = True
        else:
            materials_imported = False

        # If the element uses nodal material properties, transfer them to the nodes
        if self.element_has_nodal_properties:
            self._SetNodalProperties()

        return materials_imported

    def _SetNodalProperties(self):
        # Get fluid 1 and 2 properties
        properties_1 = self.main_model_part.Properties[1]
        properties_2 = self.main_model_part.Properties[2]

        rho_1 = properties_1.GetValue(KratosMultiphysics.DENSITY)
        rho_2 = properties_2.GetValue(KratosMultiphysics.DENSITY)
        mu_1 = properties_1.GetValue(KratosMultiphysics.DYNAMIC_VISCOSITY)
        mu_2 = properties_2.GetValue(KratosMultiphysics.DYNAMIC_VISCOSITY)

        # Check fluid 1 and 2 properties
        if rho_1 <= 0.0:
            raise Exception("DENSITY set to {0} in Properties {1}, positive number expected.".format(rho_1, properties_1.Id))
        if rho_2 <= 0.0:
            raise Exception("DENSITY set to {0} in Properties {1}, positive number expected.".format(rho_2, properties_2.Id))
        if mu_1 <= 0.0:
            raise Exception("DYNAMIC_VISCOSITY set to {0} in Properties {1}, positive number expected.".format(mu_1, properties_1.Id))
        if mu_2 <= 0.0:
            raise Exception("DYNAMIC_VISCOSITY set to {0} in Properties {1}, positive number expected.".format(mu_2, properties_2.Id))

        # Transfer density and (dynamic) viscostity to the nodes
        for node in self.main_model_part.Nodes:
            if node.GetSolutionStepValue(KratosMultiphysics.DISTANCE) <= 0.0:
                node.SetSolutionStepValue(KratosMultiphysics.DENSITY, rho_1)
                node.SetSolutionStepValue(KratosMultiphysics.DYNAMIC_VISCOSITY, mu_1)
            else:
                node.SetSolutionStepValue(KratosMultiphysics.DENSITY, rho_2)
                node.SetSolutionStepValue(KratosMultiphysics.DYNAMIC_VISCOSITY, mu_2)

    def _set_distance_function(self):
        ## Set the nodal distance function
        if (self.settings["distance_reading_settings"]["import_mode"].GetString() == "from_GiD_file"):
            DistanceUtility = DistanceImportUtility(self.main_model_part, self.settings["distance_reading_settings"])
            DistanceUtility.ImportDistance()
        elif (self.settings["distance_reading_settings"]["import_mode"].GetString() == "from_mdpa"):
            KratosMultiphysics.Logger.PrintInfo("Navier Stokes Embedded Solver","Distance function taken from the .mdpa input file.")

    def _set_level_set_convection_process(self):
        # Construct the level set convection process
        if self._bfecc_convection:
            if have_conv_diff:
                if self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 2:
                    locator = KratosMultiphysics.BinBasedFastPointLocator2D(self.main_model_part)
                    locator.UpdateSearchDatabase()
                    level_set_convection_process = KratosConvDiff.BFECCConvection2D(locator)
                else:
                    locator = KratosMultiphysics.BinBasedFastPointLocator3D(self.main_model_part)
                    locator.UpdateSearchDatabase()
                    level_set_convection_process = KratosConvDiff.BFECCConvection3D(locator)

            else:
                raise Exception("The BFECC level set convection requires the Kratos ConvectionDiffusionApplication compilation.")
        else:
            if self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 2:
                level_set_convection_process = KratosMultiphysics.LevelSetConvectionProcess2D(
                    KratosMultiphysics.DISTANCE,
                    self.main_model_part,
                    self.linear_solver)
            else:
                level_set_convection_process = KratosMultiphysics.LevelSetConvectionProcess3D(
                    KratosMultiphysics.DISTANCE,
                    KratosCFD.CONVECTIVE_VELOCITY,
                    self.main_model_part,
                    self.linear_solver,
                    0.5,    #dt_factor = 1.0
                    1.0,    #max_cfl = 1.0
                    0.7,    #cross_wind_stabilization_factor = 0.7
                    0)      #max_substeps = 0: diabled

        return level_set_convection_process

    def _set_variational_distance_process(self):
        # Construct the variational distance calculation process
        maximum_iterations = 2 #TODO: Make this user-definable
        if self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 2:
            variational_distance_process = KratosMultiphysics.VariationalDistanceCalculationProcess2D(
                self.main_model_part,
                self.linear_solver,
                maximum_iterations,
                KratosMultiphysics.VariationalDistanceCalculationProcess2D.CALCULATE_EXACT_DISTANCES_TO_PLANE)
        else:
            variational_distance_process = KratosMultiphysics.VariationalDistanceCalculationProcess3D(
                self.main_model_part,
                self.linear_solver,
                maximum_iterations,
                KratosMultiphysics.VariationalDistanceCalculationProcess3D.CALCULATE_EXACT_DISTANCES_TO_PLANE)

        return variational_distance_process

    def _set_parallel_distance_process(self):
        # Construct the parallel distance process
        if self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 2:
                locator = KratosMultiphysics.BinBasedFastPointLocator2D(self.main_model_part)
                locator.UpdateSearchDatabase()
                parallel_distance_process = KratosMultiphysics.ParallelDistanceCalculator2D()
        else:
                locator = KratosMultiphysics.BinBasedFastPointLocator3D(self.main_model_part)
                locator.UpdateSearchDatabase()
                parallel_distance_process = KratosMultiphysics.ParallelDistanceCalculator3D()

        return parallel_distance_process

    #def _set_lumped_eikonal_distance_calculation(self):
        # Construct the process for redistancing
    #    lumped_eikonal_distance_calculation = KratosCFD.LumpedEikonalDistanceCalculation(
    #            self.main_model_part,
    #            10,
    #            1.0e-8,
    #            0.1)
    #
    #    return lumped_eikonal_distance_calculation

    def _set_variational_non_eikonal_distance(self):
        #Distance re-initialization by solving a diffusion problem
        variational_non_eikonal_distance = KratosCFD.VariationalNonEikonalDistance(
                self.main_model_part,
                self.linear_solver)

        return variational_non_eikonal_distance

    def _set_distance_gradient_process(self):
        #Calculate DISTANCE_GRADIENT at nodes using ComputeNodalGradientProcess
        distance_gradient_process = KratosMultiphysics.ComputeNodalGradientProcess(
                self.main_model_part,
                KratosCFD.DISTANCE_AUX,  #KratosMultiphysics.DISTANCE,
                KratosMultiphysics.DISTANCE_GRADIENT, #KratosCFD.DISTANCE_GRADIENT_AUX,
                KratosMultiphysics.NODAL_AREA)

        return distance_gradient_process

    def _set_curvature_calculation_process(self):
        #Calculate curvature as divergence of normalized DISTANCE_GRADIENT at nodes using ComputeNodalNormalDivergenceProcess
        curvature_calculation_process = KratosMultiphysics.ComputeNodalNormalDivergenceProcess(
                self.main_model_part,
                KratosMultiphysics.DISTANCE_GRADIENT, #KratosCFD.DISTANCE_GRADIENT_AUX,
                KratosCFD.CURVATURE,
                KratosMultiphysics.NODAL_AREA)

        return curvature_calculation_process

    def _set_interface_curvature_calculation(self):
        #Calculate curvature as divergence of normalized DISTANCE_GRADIENT at nodes using LumpedInterfaceCurvatureCalculation
        interface_curvature_calculation = KratosCFD.LumpedInterfaceCurvatureCalculation(
                self.main_model_part)

        return interface_curvature_calculation

    def _set_surface_smoothing_process(self):
        #Smoothing the surface (zero DISTANCE) by solving a diffusion problem
        surface_smoothing_process = KratosCFD.SurfaceSmoothingProcess(
                self.main_model_part,
                self.linear_solver)

        return surface_smoothing_process

    def _set_mass_conservation_correction(self):
        #Restoring the distance function to compensate for the mass loss, proportional to curvature 
        #(could be velocity instead) but this is associated with the way diffusive surface smoothing process works!
        mass_conservation_correction = KratosCFD.MassConservationCorrection(
                self.main_model_part,
                True,
                "mass_conservation.log")

        return mass_conservation_correction

    def _set_interface_pressure_gradient_calculation(self):
        #Calculate pressure gradient on positive and negative sides (separately) using LumpedInterfacePositiveNegativePressureGradient
        interface_pressure_gradient_calculation = KratosCFD.LumpedInterfacePositiveNegativePressureGradient(
                self.main_model_part)

        return interface_pressure_gradient_calculation

    def _set_distance_modification_process(self):
        parameters = KratosMultiphysics.Parameters( """
        {
            "distance_threshold"                     : 1.0e-6,
            "continuous_distance"                    : true,
            "check_at_each_time_step"                : true,
            "avoid_almost_empty_elements"            : false,
            "deactivate_full_negative_elements"      : false,
            "recover_original_distance_at_each_step" : false
        }  """ )

        distance_modification_process = KratosCFD.DistanceModificationProcess(self.main_model_part, parameters)

        KratosMultiphysics.Logger.PrintInfo("DistanceModificationProcess","Construction finished.")
        distance_modification_process.ExecuteInitialize()
        distance_modification_process.ExecuteBeforeSolutionLoop()
        distance_modification_process.ExecuteInitializeSolutionStep()
        distance_modification_process.ExecuteFinalizeSolutionStep()

        return distance_modification_process

    def _set_mass_conservation_check_process(self):
        parameters = KratosMultiphysics.Parameters( """
        {
            "model_part_name"                        : "Parts_Fluid",
            "perform_corrections"                    : true,
            "correction_frequency_in_time_steps"     : 1,
            "write_to_log_file"                      : true,
            "log_file_name"                          : "mass_conservation.log"
		}  """ )

        mass_conservation_check_process = KratosCFD.MassConservationCheckProcess(self.main_model_part, parameters)

        KratosMultiphysics.Logger.PrintInfo("ApplyMassConservationCheckProcess","Construction finished.")

        first_lines_string = mass_conservation_check_process.Initialize()
        # writing first line in file
        with open("mass_conservation.log", "w") as logFile:
            logFile.write( first_lines_string )

        return mass_conservation_check_process

    def _set_find_neighbouring_elements_process(self):
        dimensions = self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
        avg_num_elements = 10
        find_neighbouring_elements_process = KratosMultiphysics.FindElementalNeighboursProcess(
                self.main_model_part, dimensions, avg_num_elements)

        return find_neighbouring_elements_process

    def _set_find_neighbouring_nodes_process(self):
        kratos_comm  = KratosMultiphysics.DataCommunicator.GetDefault()
        find_neighbouring_nodes_process = KratosMultiphysics.FindGlobalNodalNeighboursProcess(
                kratos_comm, self.main_model_part)

        return find_neighbouring_nodes_process

    def _set_hyperbolic_distance_reinitialization(self):
        distance_gradient_process_redistance = KratosMultiphysics.ComputeNodalGradientProcess(
                self.main_model_part,
                KratosMultiphysics.DISTANCE,
                KratosMultiphysics.DISTANCE_GRADIENT,
                KratosMultiphysics.NODAL_AREA)

        # Construct the process
        if self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 2:
            locator_redistance = KratosMultiphysics.BinBasedFastPointLocator2D(self.main_model_part)
            locator_redistance.UpdateSearchDatabase()
            hyperbolic_distance_reinitialization = KratosCFD.HyperbolicDistanceReinitialization2D(
                self.main_model_part, locator_redistance, distance_gradient_process_redistance,
                50000, 1.0e-9, 1.0e-6)
        else:
            locator_redistance = KratosMultiphysics.BinBasedFastPointLocator3D(self.main_model_part)
            locator_redistance.UpdateSearchDatabase()
            hyperbolic_distance_reinitialization = KratosCFD.HyperbolicDistanceReinitialization3D(
                self.main_model_part, locator_redistance, distance_gradient_process_redistance,
                50000, 1.0e-9, 1.0e-6)

        return hyperbolic_distance_reinitialization


    def _set_contact_angle_evaluator(self):
        contact_angle_evaluator = KratosCFD.ContactAngleEvaluator(self.main_model_part)

        return contact_angle_evaluator