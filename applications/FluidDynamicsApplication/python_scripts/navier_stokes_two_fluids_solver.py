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

from pathlib import Path

def CreateSolver(model, custom_settings):
    return NavierStokesTwoFluidsSolver(model, custom_settings)

class NavierStokesTwoFluidsSolver(FluidSolver):

    @classmethod
    def GetDefaultParameters(cls):
        ##settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
        {
            "solver_type": "two_fluids",
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
            "time_scheme": "bdf2",
            "compute_reactions": false,
            "analysis_type": "non_linear",
            "reform_dofs_at_each_step": false,
            "consider_periodic_conditions": false,
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
            "acceleration_limitation": true,
            "formulation": {
                "dynamic_tau": 1.0,
                "surface_tension": false,
                "mass_source":false,
                "momentum_correction":false
            },
            "levelset_convection_settings": {
                "max_CFL" : 1.0,
                "max_substeps" : 0,
                "eulerian_error_compensation" : false,
                "element_type" : "levelset_convection_supg",
                "element_settings" : {
                    "dynamic_tau" : 0.0,
                    "cross_wind_stabilization_factor" : 0.7
                }
            },
            "distance_reinitialization": "variational",
            "parallel_redistance_max_layers" : 25,
            "distance_smoothing": false,
            "distance_smoothing_coefficient": 1.0,
            "distance_modification_settings": {
                "model_part_name": "",
                "distance_threshold": 1e-5,
                "continuous_distance": true,
                "check_at_each_time_step": true,
                "avoid_almost_empty_elements": false,
                "deactivate_full_negative_elements": false
            }
        }""")

        default_settings.AddMissingParameters(super(NavierStokesTwoFluidsSolver, cls).GetDefaultParameters())
        return default_settings

    def __init__(self, model, custom_settings):
        # TODO: DO SOMETHING IN HERE TO REMOVE THE "time_order" FROM THE DEFAULT SETTINGS BUT KEEPING THE BACKWARDS COMPATIBILITY

        if custom_settings.Has("levelset_convection_settings"):
            if custom_settings["levelset_convection_settings"].Has("levelset_splitting"):
                custom_settings["levelset_convection_settings"].RemoveValue("levelset_splitting")
                KratosMultiphysics.Logger.PrintWarning("NavierStokesTwoFluidsSolver", "\'levelset_splitting\' has been temporarily deactivated. Using the standard levelset convection with no splitting.")

        #TODO: Remove this after the retrocompatibility period
        if custom_settings.Has("bfecc_convection"):
            KratosMultiphysics.Logger.PrintWarning("NavierStokesTwoFluidsSolver", "the semi-Lagrangian \'bfecc_convection\' is no longer supported. Using the standard Eulerian levelset convection.")
            custom_settings.RemoveValue("bfecc_convection")
            if custom_settings.Has("bfecc_number_substeps"):
                custom_settings.RemoveValue("bfecc_number_substeps")

        super(NavierStokesTwoFluidsSolver,self).__init__(model,custom_settings)

        self.element_name = "TwoFluidNavierStokes"
        self.condition_name = "TwoFluidNavierStokesWallCondition"
        self.element_integrates_in_time = True
        self.element_has_nodal_properties = True

        self.min_buffer_size = 3

        # Set the levelset characteristic variables and add them to the convection settings
        # These are required to be set as some of the auxiliary processes admit user-defined variables
        self._levelset_variable = KratosMultiphysics.DISTANCE
        self._levelset_gradient_variable = KratosMultiphysics.DISTANCE_GRADIENT
        self._levelset_convection_variable = KratosMultiphysics.VELOCITY
        self.settings["levelset_convection_settings"].AddEmptyValue("levelset_variable_name").SetString("DISTANCE")
        self.settings["levelset_convection_settings"].AddEmptyValue("levelset_gradient_variable_name").SetString("DISTANCE_GRADIENT")
        self.settings["levelset_convection_settings"].AddEmptyValue("levelset_convection_variable_name").SetString("VELOCITY")

        dynamic_tau = self.settings["formulation"]["dynamic_tau"].GetDouble()
        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DYNAMIC_TAU, dynamic_tau)

        surface_tension = False
        if (self.settings["formulation"].Has("surface_tension")):
            surface_tension = self.settings["formulation"]["surface_tension"].GetBool()
        self.main_model_part.ProcessInfo.SetValue(KratosCFD.SURFACE_TENSION, surface_tension)

        self.momentum_correction = False
        if self.settings["formulation"].Has("momentum_correction"):
            self.momentum_correction = self.settings["formulation"]["momentum_correction"].GetBool()
        self.main_model_part.ProcessInfo.SetValue(KratosCFD.MOMENTUM_CORRECTION, self.momentum_correction)

        self._reinitialization_type = self.settings["distance_reinitialization"].GetString()

        self._distance_smoothing = self.settings["distance_smoothing"].GetBool()
        smoothing_coefficient = self.settings["distance_smoothing_coefficient"].GetDouble()
        self.main_model_part.ProcessInfo.SetValue(KratosCFD.SMOOTHING_COEFFICIENT, smoothing_coefficient)

        self._apply_acceleration_limitation = self.settings["acceleration_limitation"].GetBool()

        ## Set the distance reading filename
        # TODO: remove the manual "distance_file_name" set as soon as the problem type one has been tested.
        if (self.settings["distance_reading_settings"]["import_mode"].GetString() == "from_GiD_file"):
            self.settings["distance_reading_settings"]["distance_file_name"].SetString(self.settings["model_import_settings"]["input_filename"].GetString()+".post.res")

        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Construction of NavierStokesTwoFluidsSolver finished.")

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
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE_GRADIENT)     # Distance gradient nodal values

        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Fluid solver variables added correctly.")

    def PrepareModelPart(self):
        # Initialize the level-set function
        if not self.main_model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED]:
            ## Setting the nodal distance
            self.__SetDistanceFunction()

        # Call the base solver PrepareModelPart()
        super(NavierStokesTwoFluidsSolver, self).PrepareModelPart()

    def Initialize(self):
        computing_model_part = self.GetComputingModelPart()
        # Calculate boundary normals
        KratosMultiphysics.NormalCalculationUtils().CalculateOnSimplex(
            computing_model_part,
            computing_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE])

        # Finding nodal and elemental neighbors
        data_communicator = computing_model_part.GetCommunicator().GetDataCommunicator()
        neighbour_search = KratosMultiphysics.FindGlobalNodalNeighboursProcess(
            data_communicator,
            computing_model_part)
        neighbour_search.Execute()

        elemental_neighbour_search = KratosMultiphysics.GenericFindElementalNeighboursProcess(
            computing_model_part)
        elemental_neighbour_search.Execute()

        # Set and initialize the solution strategy
        solution_strategy = self._GetSolutionStrategy()
        solution_strategy.SetEchoLevel(self.settings["echo_level"].GetInt())
        solution_strategy.Initialize()

        # Set nodal properties after setting distance(level-set).
        self._SetNodalProperties()

        # Initialize the distance correction process
        self._GetDistanceModificationProcess().ExecuteInitialize()
        self._GetDistanceModificationProcess().ExecuteInitializeSolutionStep()

        #Here the initial water volume of the system is calculated without considering inlet and outlet flow rate
        self.initial_system_volume=KratosCFD.FluidAuxiliaryUtilities.CalculateFluidNegativeVolume(self.GetComputingModelPart())

        # Instantiate the level set convection process
        # Note that is is required to do this in here in order to validate the defaults and set the corresponding distance gradient flag
        # Note that the nodal gradient of the distance is required either for the eulerian BFECC limiter or by the algebraic element antidiffusivity
        self._GetLevelSetConvectionProcess()

        self.mass_source = False
        if self.settings["formulation"].Has("mass_source"):
            self.mass_source = self.settings["formulation"]["mass_source"].GetBool()

        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Solver initialization finished.")

    def InitializeSolutionStep(self):

        if self.momentum_correction:
            KratosMultiphysics.VariableUtils().SetNonHistoricalVariable(KratosCFD.DISTANCE_CORRECTION, 0.0, self.main_model_part.Nodes)

        # Inlet and outlet water discharge is calculated for current time step, first discharge and the considering the time step inlet and outlet volume is calculated
        if self.mass_source:
            outlet_discharge = KratosCFD.FluidAuxiliaryUtilities.CalculateFlowRateNegativeSkin(self.GetComputingModelPart(),KratosMultiphysics.OUTLET)
            inlet_discharge = KratosCFD.FluidAuxiliaryUtilities.CalculateFlowRateNegativeSkin(self.GetComputingModelPart(),KratosMultiphysics.INLET)
            current_dt = self.main_model_part.ProcessInfo[KratosMultiphysics.DELTA_TIME]
            inlet_volume = -current_dt * inlet_discharge
            outlet_volume = current_dt * outlet_discharge

            # System water volume is calculated for current time step considering inlet and outlet discharge.
            system_volume = inlet_volume + self.initial_system_volume - outlet_volume

        # Recompute the BDF2 coefficients
        (self.time_discretization).ComputeAndSaveBDFCoefficients(self.GetComputingModelPart().ProcessInfo)

        # Perform the level-set convection according to the previous step velocity
        self._PerformLevelSetConvection()

        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Level-set convection is performed.")

        # filtering noises is necessary for curvature calculation
        if (self._distance_smoothing):
            # distance gradient is used as a boundary condition for smoothing process
            self._GetDistanceGradientProcess().Execute()
            self._GetDistanceSmoothingProcess().Execute()
            KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Smoothing process is finished.")

        if (self.main_model_part.ProcessInfo[KratosCFD.SURFACE_TENSION]):
            # distance gradient is called again to comply with the smoothed/modified DISTANCE
            self._GetDistanceGradientProcess().Execute()
            # curvature is calculated using nodal distance gradient
            self._GetDistanceCurvatureProcess().Execute()
            # it is needed to store level-set consistent nodal PRESSURE_GRADIENT for stabilization purpose
            self._GetConsistentNodalPressureGradientProcess().Execute()

        # TODO: Performing mass conservation check and correction process

        # Perform distance correction to prevent ill-conditioned cuts
        self._GetDistanceModificationProcess().ExecuteInitializeSolutionStep()

        # Update the DENSITY and DYNAMIC_VISCOSITY values according to the new level-set
        self._SetNodalProperties()

        # Initialize the solver current step
        self._GetSolutionStrategy().InitializeSolutionStep()

        # Accumulative water volume error ratio due to level set. Adding source term
        if self.mass_source:
            water_volume_after_transport = KratosCFD.FluidAuxiliaryUtilities.CalculateFluidNegativeVolume(self.GetComputingModelPart())
            volume_error = (water_volume_after_transport - system_volume) / system_volume
            self.initial_system_volume=system_volume
        else:
            volume_error=0

        self.main_model_part.ProcessInfo.SetValue(KratosCFD.VOLUME_ERROR, volume_error)

        # We set this value at every time step as other processes/solvers also use them
        dynamic_tau = self.settings["formulation"]["dynamic_tau"].GetDouble()
        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DYNAMIC_TAU, dynamic_tau)


    def FinalizeSolutionStep(self):
        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Mass and momentum conservation equations are solved.")

        # Recompute the distance field according to the new level-set position
        if self._reinitialization_type != "none":
            self._GetDistanceReinitializationProcess().Execute()
            KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Redistancing process is finished.")

        # Prepare distance correction for next step
        self._GetDistanceModificationProcess().ExecuteFinalizeSolutionStep()

        # Finalize the solver current step
        self._GetSolutionStrategy().FinalizeSolutionStep()
        # Limit the obtained acceleration for the next step
        # This limitation should be called on the second solution step onwards (e.g. STEP=2 for BDF2)
        # We intentionally avoid correcting the acceleration in the first resolution step as this might cause problems with zero initial conditions
        if self._apply_acceleration_limitation and self.main_model_part.ProcessInfo[KratosMultiphysics.STEP] > 1:
            self._GetAccelerationLimitationUtility().Execute()

    def _PerformLevelSetConvection(self):
        # Solve the levelset convection problem
        self._GetLevelSetConvectionProcess().Execute()

    # TODO: Remove this method as soon as the subproperties are available
    def _SetPhysicalProperties(self):
        warn_msg  = '\nThe materials import mechanism used in the two fluids solver is DEPRECATED!\n'
        warn_msg += 'It will be removed to use the base fluid_solver.py one as soon as the subproperties are available.\n'
        KratosMultiphysics.Logger.PrintWarning('\n\x1b[1;31mDEPRECATION-WARNING\x1b[0m', warn_msg)

        # Check if the fluid properties are provided using a .json file
        materials_filename = self.settings["material_import_settings"]["materials_filename"].GetString()
        if (materials_filename != ""):
            data_comm = KratosMultiphysics.ParallelEnvironment.GetDefaultDataCommunicator() # only using the global comm as the Communicators are not yet created when running in MPI. Hotfix since this method will disappear completely when using subproperties!

            def GetAuxMaterialsFileName(mat_file_name, prop_id):
                p_mat_file_name = Path(mat_file_name)
                new_stem = "{}_p{}".format(p_mat_file_name.stem, prop_id)
                return str(p_mat_file_name.with_name(new_stem).with_suffix(p_mat_file_name.suffix))

            with open(materials_filename,'r') as materials_file:
                materials = KratosMultiphysics.Parameters(materials_file.read())

            if data_comm.Rank() == 0:
                # Create and read an auxiliary materials file for each one of the fields (only on one rank)
                for i_material in materials["properties"].values():
                    aux_materials = KratosMultiphysics.Parameters()
                    aux_materials.AddEmptyArray("properties")
                    aux_materials["properties"].Append(i_material)

                    aux_materials_filename = GetAuxMaterialsFileName(materials_filename, i_material["properties_id"].GetInt())
                    with open(aux_materials_filename,'w') as aux_materials_file:
                        aux_materials_file.write(aux_materials.WriteJsonString())

            data_comm.Barrier()

            # read the files on all ranks
            for i_material in materials["properties"].values():
                aux_materials_filename = GetAuxMaterialsFileName(materials_filename, i_material["properties_id"].GetInt())
                aux_material_settings = KratosMultiphysics.Parameters("""{"Parameters": {"materials_filename": ""}} """)
                aux_material_settings["Parameters"]["materials_filename"].SetString(aux_materials_filename)
                KratosMultiphysics.ReadMaterialsUtility(aux_material_settings, self.model)

            data_comm.Barrier()

            if data_comm.Rank() == 0:
                # remove aux files after every rank read them
                for i_material in materials["properties"].values():
                    aux_materials_filename = GetAuxMaterialsFileName(materials_filename, i_material["properties_id"].GetInt())
                    KratosUtilities.DeleteFileIfExisting(aux_materials_filename)

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
            if node.GetSolutionStepValue(self._levelset_variable) <= 0.0:
                node.SetSolutionStepValue(KratosMultiphysics.DENSITY, rho_1)
                node.SetSolutionStepValue(KratosMultiphysics.DYNAMIC_VISCOSITY, mu_1)
            else:
                node.SetSolutionStepValue(KratosMultiphysics.DENSITY, rho_2)
                node.SetSolutionStepValue(KratosMultiphysics.DYNAMIC_VISCOSITY, mu_2)

    def __SetDistanceFunction(self):
        ## Set the nodal distance function
        if (self.settings["distance_reading_settings"]["import_mode"].GetString() == "from_GiD_file"):
            DistanceUtility = DistanceImportUtility(self.main_model_part, self.settings["distance_reading_settings"])
            DistanceUtility.ImportDistance()
        elif (self.settings["distance_reading_settings"]["import_mode"].GetString() == "from_mdpa"):
            KratosMultiphysics.Logger.PrintInfo("Navier Stokes Embedded Solver","Distance function taken from the .mdpa input file.")

    def _GetAccelerationLimitationUtility(self):
        if not hasattr(self, '_acceleration_limitation_utility'):
            self._acceleration_limitation_utility = self.__CreateAccelerationLimitationUtility()
        return self._acceleration_limitation_utility

    def _GetRedistancingLinearSolver(self):
        # A linear solver configured specifically for distance re-initialization process
        if not hasattr(self, '_redistancing_linear_solver'):
            self._redistancing_linear_solver = self._CreateLinearSolver() # TODO: add customized configuration
        return self._redistancing_linear_solver

    def _GetLevelsetLinearSolver(self):
        # A linear solver configured specifically for the level-set convection process
        if not hasattr(self, '_levelset_linear_solver'):
            self._levelset_linear_solver = self._CreateLinearSolver() # TODO: add customized configuration
        return self._levelset_linear_solver

    def _GetSmoothingLinearSolver(self):
        # A linear solver configured specifically for the distance smoothing process
        if not hasattr(self, '_smoothing_linear_solver'):
            self._smoothing_linear_solver = self._CreateLinearSolver() # TODO: add customized configuration
        return self._smoothing_linear_solver

    def _GetLevelSetConvectionProcess(self):
        if not hasattr(self, '_level_set_convection_process'):
            self._level_set_convection_process = self._CreateLevelSetConvectionProcess()
        return self._level_set_convection_process

    def _GetDistanceReinitializationProcess(self):
        if not hasattr(self, '_distance_reinitialization_process'):
            self._distance_reinitialization_process = self._CreateDistanceReinitializationProcess()
        return self._distance_reinitialization_process

    def _GetDistanceSmoothingProcess(self):
        if not hasattr(self, '_distance_smoothing_process'):
            self._distance_smoothing_process = self._CreateDistanceSmoothingProcess()
        return self._distance_smoothing_process

    def _GetDistanceGradientProcess(self):
        if not hasattr(self, '_distance_gradient_process'):
            self._distance_gradient_process = self._CreateDistanceGradientProcess()
        return self._distance_gradient_process

    def _GetDistanceCurvatureProcess(self):
        if not hasattr(self, '_distance_curvature_process'):
            self._distance_curvature_process = self._CreateDistanceCurvatureProcess()
        return self._distance_curvature_process

    def _GetConsistentNodalPressureGradientProcess(self):
        if not hasattr(self, '_consistent_nodal_pressure_gradient_process'):
            self._consistent_nodal_pressure_gradient_process = self._CreateConsistentNodalPressureGradientProcess()
        return self._consistent_nodal_pressure_gradient_process

    def _GetDistanceModificationProcess(self):
        if not hasattr(self, '_distance_modification_process'):
            self._distance_modification_process = self.__CreateDistanceModificationProcess()
        return self._distance_modification_process

    def __CreateAccelerationLimitationUtility(self):
        maximum_multiple_of_g_acceleration_allowed = 5.0
        acceleration_limitation_utility = KratosCFD.AccelerationLimitationUtilities(
            self.GetComputingModelPart(),
            maximum_multiple_of_g_acceleration_allowed)

        return acceleration_limitation_utility

    def _CreateLevelSetConvectionProcess(self):
        # Construct the level set convection process
        domain_size = self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
        computing_model_part = self.GetComputingModelPart()
        linear_solver = self._GetLevelsetLinearSolver()
        levelset_convection_settings = self.settings["levelset_convection_settings"]
        if domain_size == 2:
            level_set_convection_process = KratosMultiphysics.LevelSetConvectionProcess2D(
                computing_model_part,
                linear_solver,
                levelset_convection_settings)
        else:
            level_set_convection_process = KratosMultiphysics.LevelSetConvectionProcess3D(
                computing_model_part,
                linear_solver,
                levelset_convection_settings)

        return level_set_convection_process

    def _CreateDistanceReinitializationProcess(self):
        # Construct the variational distance calculation process
        if (self._reinitialization_type == "variational"):
            maximum_iterations = 2 #TODO: Make this user-definable
            linear_solver = self._GetRedistancingLinearSolver()
            computing_model_part = self.GetComputingModelPart()
            if self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 2:
                distance_reinitialization_process = KratosMultiphysics.VariationalDistanceCalculationProcess2D(
                    computing_model_part,
                    linear_solver,
                    maximum_iterations,
                    KratosMultiphysics.VariationalDistanceCalculationProcess2D.CALCULATE_EXACT_DISTANCES_TO_PLANE)
            else:
                distance_reinitialization_process = KratosMultiphysics.VariationalDistanceCalculationProcess3D(
                    computing_model_part,
                    linear_solver,
                    maximum_iterations,
                    KratosMultiphysics.VariationalDistanceCalculationProcess3D.CALCULATE_EXACT_DISTANCES_TO_PLANE)

        elif (self._reinitialization_type == "parallel"):
            #TODO: move all this to solver settings
            layers = self.settings["parallel_redistance_max_layers"].GetInt()
            parallel_distance_settings = KratosMultiphysics.Parameters("""{
                "max_levels" : 25,
                "max_distance" : 1.0,
                "calculate_exact_distances_to_plane" : true
            }""")
            parallel_distance_settings["max_levels"].SetInt(layers)
            if self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 2:
                distance_reinitialization_process = KratosMultiphysics.ParallelDistanceCalculationProcess2D(
                    self.main_model_part,
                    parallel_distance_settings)
            else:
                distance_reinitialization_process = KratosMultiphysics.ParallelDistanceCalculationProcess3D(
                    self.main_model_part,
                    parallel_distance_settings)
        elif (self._reinitialization_type == "none"):
                KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Redistancing is turned off.")
        else:
            raise Exception("Please use a valid distance reinitialization type or set it as \'none\'. Valid types are: \'variational\' and \'parallel\'.")

        return distance_reinitialization_process

    def _CreateDistanceSmoothingProcess(self):
        # construct the distance smoothing process
        linear_solver = self._GetSmoothingLinearSolver()
        if self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 2:
            distance_smoothing_process = KratosCFD.DistanceSmoothingProcess2D(
            self.main_model_part,
            linear_solver)
        else:
            distance_smoothing_process = KratosCFD.DistanceSmoothingProcess3D(
            self.main_model_part,
            linear_solver)

        return distance_smoothing_process

    def _CreateDistanceGradientProcess(self):
        distance_gradient_process = KratosMultiphysics.ComputeNodalGradientProcess(
                self.main_model_part,
                self._levelset_variable,
                self._levelset_gradient_variable,
                KratosMultiphysics.NODAL_AREA)

        return distance_gradient_process

    def _CreateDistanceCurvatureProcess(self):
        distance_curvature_process = KratosMultiphysics.ComputeNonHistoricalNodalNormalDivergenceProcess(
                self.main_model_part,
                self._levelset_gradient_variable,
                KratosCFD.CURVATURE,
                KratosMultiphysics.NODAL_AREA)

        return distance_curvature_process

    def _CreateConsistentNodalPressureGradientProcess(self):
        consistent_nodal_pressure_gradient_process = KratosCFD.CalulateLevelsetConsistentNodalGradientProcess(
                self.main_model_part)

        return consistent_nodal_pressure_gradient_process

    def __CreateDistanceModificationProcess(self):
        # Set suitable distance correction settings for free-surface problems
        # Note that the distance modification process is applied to the computing model part
        distance_modification_settings = self.settings["distance_modification_settings"]
        distance_modification_settings.ValidateAndAssignDefaults(self.GetDefaultParameters()["distance_modification_settings"])
        distance_modification_settings["model_part_name"].SetString(self.GetComputingModelPart().FullName())

        # Check user provided settings
        if not distance_modification_settings["continuous_distance"].GetBool():
            distance_modification_settings["continuous_distance"].SetBool(True)
            KratosMultiphysics.Logger.PrintWarning("Provided distance correction \'continuous_distance\' is \'False\'. Setting to \'True\'.")
        if not distance_modification_settings["check_at_each_time_step"].GetBool():
            distance_modification_settings["check_at_each_time_step"].SetBool(True)
            KratosMultiphysics.Logger.PrintWarning("Provided distance correction \'check_at_each_time_step\' is \'False\'. Setting to \'True\'.")
        if distance_modification_settings["avoid_almost_empty_elements"].GetBool():
            distance_modification_settings["avoid_almost_empty_elements"].SetBool(False)
            KratosMultiphysics.Logger.PrintWarning("Provided distance correction \'avoid_almost_empty_elements\' is \'True\'. Setting to \'False\' to avoid modifying the distance sign.")
        if distance_modification_settings["deactivate_full_negative_elements"].GetBool():
            distance_modification_settings["deactivate_full_negative_elements"].SetBool(False)
            KratosMultiphysics.Logger.PrintWarning("Provided distance correction \'deactivate_full_negative_elements\' is \'True\'. Setting to \'False\' to avoid deactivating the negative volume (e.g. water).")

        # Create and return the distance correction process
        return KratosCFD.DistanceModificationProcess(
            self.model,
            distance_modification_settings)
