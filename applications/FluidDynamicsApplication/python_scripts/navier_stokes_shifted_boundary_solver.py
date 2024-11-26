# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD

# Import base class file
from KratosMultiphysics.FluidDynamicsApplication.fluid_solver import FluidSolver


class ShiftedBoundaryFormulation(object):
    """Helper class to define shifted boundary dependent parameters."""
    def __init__(self, formulation_settings):
        self.element_name = None
        self.condition_name = None
        self.process_info_data = {}
        self.element_has_nodal_properties = False
        self.historical_nodal_properties_variables_list = []
        self.non_historical_nodal_properties_variables_list = []

        if formulation_settings.Has("element_type"):
            element_type = formulation_settings["element_type"].GetString()
            if element_type == "shifted_boundary_weakly_compressible_navier_stokes":
                self._SetUpShiftedBoundaryWeaklyCompressibleNavierStokes(formulation_settings)
        else:
            raise RuntimeError("Argument \'element_type\' not found in formulation settings.")

    def SetProcessInfo(self, model_part):
        for variable,value in self.process_info_data.items():
            model_part.ProcessInfo[variable] = value

    def _SetUpShiftedBoundaryWeaklyCompressibleNavierStokes(self, formulation_settings):
        default_settings = KratosMultiphysics.Parameters(r"""{
            "element_type": "shifted_boundary_weakly_compressible_navier_stokes",
            "is_slip": false,
            "slip_length": 1.0e8,
            "penalty_coefficient": 10.0,
            "dynamic_tau": 1.0,
            "level_set_type": "point-based",
            "conforming_basis" : true,
            "extension_operator_type" : "MLS",
            "mls_extension_operator_order" : 1
        }""")
        formulation_settings.ValidateAndAssignDefaults(default_settings)

        self.element_name = "ShiftedBoundaryWeaklyCompressibleNavierStokes"
        self.condition_name = "NavierStokesWallCondition"
        self.sbm_interface_condition_name = "ShiftedBoundaryWallCondition"
        self.level_set_type = formulation_settings["level_set_type"].GetString()
        # Error that discontinuous is not supported yet
        if self.level_set_type != "continuous" and self.level_set_type != "discontinuous" and self.level_set_type != "point-based":
            err_msg = 'Provided level set type is unknown. Available types for MLS-based SBM are \'continuous\', \'discontinuous\' and \'point-based\'.'
            raise Exception(err_msg)
        self.element_integrates_in_time = True
        self.element_has_nodal_properties = True
        self.historical_nodal_properties_variables_list = [KratosMultiphysics.DENSITY]
        self.non_historical_nodal_properties_variables_list = [KratosMultiphysics.SOUND_VELOCITY]

        self.process_info_data[KratosMultiphysics.DYNAMIC_TAU] = formulation_settings["dynamic_tau"].GetDouble()
        self.process_info_data[KratosMultiphysics.PENALTY_COEFFICIENT] = formulation_settings["penalty_coefficient"].GetDouble()
        if formulation_settings["is_slip"].GetBool():
            self.process_info_data[KratosCFD.SLIP_LENGTH] = formulation_settings["slip_length"].GetDouble()


def CreateSolver(model, custom_settings):
    return NavierStokesShiftedBoundaryMonolithicSolver(model, custom_settings)

class NavierStokesShiftedBoundaryMonolithicSolver(FluidSolver):

    def __GetDistanceModificationDefaultSettings(self, level_set_type):
        if level_set_type == "continuous":
            return self.__GetContinuousDistanceModificationDefaultSettings()
        elif level_set_type == "discontinuous":
            return self.__GetDiscontinuousDistanceModificationDefaultSettings()
        else:
            err_msg = 'Provided level set type is: \'' + level_set_type + '\'. Expected \'continuous\' or \'discontinuous\'.'
            raise Exception(err_msg)

    @classmethod
    def __GetContinuousDistanceModificationDefaultSettings(cls):
        return KratosMultiphysics.Parameters(r'''{
            "model_part_name": "",
            "distance_threshold": 1.0e-12,
            "continuous_distance": true,
            "check_at_each_time_step": false,
            "avoid_almost_empty_elements": false,
            "deactivate_full_negative_elements": true
        }''')

    @classmethod
    def __GetDiscontinuousDistanceModificationDefaultSettings(cls):
        return KratosMultiphysics.Parameters(r'''{
            "model_part_name": "",
            "distance_threshold": 1.0e-12,
            "continuous_distance": false,
            "check_at_each_time_step": false,
            "avoid_almost_empty_elements": false,
            "deactivate_full_negative_elements": false
        }''')

    @classmethod
    def GetDefaultParameters(cls):
        ##settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
        {
            "solver_type": "shifted_boundary_from_defaults",
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
            "distance_modification_settings": {
            },
            "maximum_iterations": 7,
            "echo_level": 0,
            "time_order": 2,
            "time_scheme": "bdf2",
            "compute_reactions": false,
            "analysis_type": "non_linear",
            "reform_dofs_at_each_step": false,
            "consider_periodic_conditions": false,
            "assign_neighbour_elements_to_conditions": true,
            "relative_velocity_tolerance": 1e-3,
            "absolute_velocity_tolerance": 1e-5,
            "relative_pressure_tolerance": 1e-3,
            "absolute_pressure_tolerance": 1e-5,
            "linear_solver_settings"       : {
                "solver_type"         : "amgcl"
            },
            "volume_model_part_name" : "volume_model_part",
            "skin_parts": [""],
            "no_skin_parts":[""],
            "time_stepping": {
                "automatic_time_step" : true,
                "CFL_number"          : 1,
                "minimum_delta_time"  : 1e-2,
                "maximum_delta_time"  : 1.0,
                "time_step"           : 0.0
            },
            "move_mesh_flag": false,
            "formulation": {
                "element_type": "shifted_boundary_weakly_compressible_navier_stokes"
            }
        }""")

        default_settings.AddMissingParameters(super(NavierStokesShiftedBoundaryMonolithicSolver, cls).GetDefaultParameters())
        return default_settings

    def __init__(self, model, custom_settings):
        super(NavierStokesShiftedBoundaryMonolithicSolver,self).__init__(model,custom_settings)

        self.min_buffer_size = 3
        self.shifted_boundary_formulation = ShiftedBoundaryFormulation(self.settings["formulation"])
        self.element_name = self.shifted_boundary_formulation.element_name
        self.condition_name = self.shifted_boundary_formulation.condition_name
        self.sbm_interface_condition_name = self.shifted_boundary_formulation.sbm_interface_condition_name + str(self.settings["domain_size"].GetInt()) + "D"
        self.level_set_type = self.shifted_boundary_formulation.level_set_type
        self.element_integrates_in_time = self.shifted_boundary_formulation.element_integrates_in_time
        self.element_has_nodal_properties = self.shifted_boundary_formulation.element_has_nodal_properties
        self.historical_nodal_properties_variables_list = self.shifted_boundary_formulation.historical_nodal_properties_variables_list
        self.non_historical_nodal_properties_variables_list = self.shifted_boundary_formulation.non_historical_nodal_properties_variables_list

        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Construction of NavierStokesShiftedBoundaryMonolithicSolver finished.")

    def AddVariables(self):
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PRESSURE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ACCELERATION)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.MESH_VELOCITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.BODY_FORCE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_AREA)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION_WATER_PRESSURE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.EXTERNAL_PRESSURE)
        if self.level_set_type == "continuous" or self.level_set_type == "discontinuous":
            self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE)              # Distance function nodal values
            self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE_GRADIENT)     # Distance gradient nodal values
        self.main_model_part.AddNodalSolutionStepVariable(KratosCFD.EMBEDDED_WET_PRESSURE)          # Post-process variable (stores the fluid nodes pressure and is set to 0 in the structure ones)
        self.main_model_part.AddNodalSolutionStepVariable(KratosCFD.EMBEDDED_WET_VELOCITY)          # Post-process variable (stores the fluid nodes velocity and is set to 0 in the structure ones)

        # Adding variables required for the nodal material properties
        if self.element_has_nodal_properties:
            for variable in self.historical_nodal_properties_variables_list:
                self.main_model_part.AddNodalSolutionStepVariable(variable)

        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Shifted-boundary fluid solver variables added correctly.")

    def AddDofs(self):
        # Add formulation DOFs and reactions
        super().AddDofs()

    def PrepareModelPart(self):
        # Call the base solver PrepareModelPart()
        super(NavierStokesShiftedBoundaryMonolithicSolver, self).PrepareModelPart()

        # Set the extra requirements of the shifted-boundary formulation
        if not self.main_model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED]:
            ## Sets the shifted-boundary formulation configuration
            self.__SetShiftedBoundaryFormulation()
            ## Setting the nodal distance
            if self.level_set_type == "continuous" or self.level_set_type == "discontinuous":
                self.__SetDistanceFunction()

    def Initialize(self):
        # If the solver requires an instance of the stabilized shifted boundary formulation class, set the process info variables
        #TODO??
        if hasattr(self, 'shifted_boundary_formulation'):
            self.shifted_boundary_formulation.SetProcessInfo(self.GetComputingModelPart())

        # Construct and initialize the solution strategy
        #TODO "Error: Constitutive Law not initialized for Element ShiftedBoundaryFluidElement #105033"
        # if strategy is initialized after set up of interface utility (deactivation of elements?)
        solution_strategy = self._GetSolutionStrategy()
        solution_strategy.SetEchoLevel(self.settings["echo_level"].GetInt())
        solution_strategy.Initialize()

        # # Set the distance modification process
        # self.GetDistanceModificationProcess().ExecuteInitialize()
        # # Correct the distance field
        # self.GetDistanceModificationProcess().ExecuteInitializeSolutionStep()
        #TODO OR
        # Avoid zeros with positive epsilon
        # TODO remove DISTANCEs
        if self.level_set_type == "continuous" or self.level_set_type == "discontinuous":
            tol = 1.0e-10
            for node in self.GetComputingModelPart().Nodes:
                dist = node.GetSolutionStepValue(KratosMultiphysics.DISTANCE)
                if abs(dist) < tol:
                    node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, 0, tol)

        # Create shifted-boundary meshless interface utility and calculate extension operator requiring nodal and elemental neighbors
        self.__SetUpInterfaceUtility()

        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Solver initialization finished.")

    def InitializeSolutionStep(self):
        # Compute the BDF coefficients
        (self.time_discretization).ComputeAndSaveBDFCoefficients(self.GetComputingModelPart().ProcessInfo)

        # Call the base solver InitializeSolutionStep()
        super(NavierStokesShiftedBoundaryMonolithicSolver, self).InitializeSolutionStep()

    def _SetNodalProperties(self):
        set_density = KratosMultiphysics.DENSITY in self.historical_nodal_properties_variables_list
        set_sound_velocity = KratosMultiphysics.SOUND_VELOCITY in self.non_historical_nodal_properties_variables_list

        # Get density and dynamic viscosity from the properties of the first element
        for ele in self.main_model_part.Elements:
            # Get DENSITY from properties
            if set_density:
                rho = ele.Properties.GetValue(KratosMultiphysics.DENSITY)
                if rho <= 0.0:
                    raise Exception("DENSITY set to {0} in Properties {1}, positive number expected.".format(rho,ele.Properties.Id))
            # Get SOUND_VELOCITY
            if set_sound_velocity:
                if ele.Properties.Has(KratosMultiphysics.SOUND_VELOCITY):
                    sound_velocity = ele.Properties.GetValue(KratosMultiphysics.SOUND_VELOCITY)
                else:
                    sound_velocity = 1.0e+12 # Default sound velocity value
                    KratosMultiphysics.Logger.PrintWarning('No \'SOUND_VELOCITY\' value found in Properties {0}. Setting default value {1}'.format(ele.Properties.Id, sound_velocity))
                if sound_velocity <= 0.0:
                    raise Exception("SOUND_VELOCITY set to {0} in Properties {1}, positive number expected.".format(sound_velocity, ele.Properties.Id))
            break
        else:
            raise Exception("No fluid elements found in the main model part.")

        # Transfer the obtained properties to the nodes
        if set_density:
            KratosMultiphysics.VariableUtils().SetVariable(KratosMultiphysics.DENSITY, rho, self.main_model_part.Nodes)
        if set_sound_velocity:
            KratosMultiphysics.VariableUtils().SetNonHistoricalVariable(KratosMultiphysics.SOUND_VELOCITY, sound_velocity, self.main_model_part.Nodes)

    def __SetShiftedBoundaryFormulation(self):
        # Set the SLIP elemental flag
        if (self.settings["formulation"]["is_slip"].GetBool()):
            KratosMultiphysics.VariableUtils().SetFlag(KratosMultiphysics.SLIP, True, self.GetComputingModelPart().Elements)
        else:
            # Set the SLIP elemental flag to false in the entire domain
            KratosMultiphysics.VariableUtils().SetFlag(KratosMultiphysics.SLIP, False, self.GetComputingModelPart().Elements)

        # Save the formulation settings in the ProcessInfo
        self.shifted_boundary_formulation.SetProcessInfo(self.main_model_part)

    def __SetDistanceFunction(self):
        ## Set the nodal distance function
        if (self.settings["distance_reading_settings"]["import_mode"].GetString() == "from_mdpa"):
            KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__,"Distance function taken from the .mdpa input file.")
            # Recall to swap the distance sign (GiD considers d<0 in the fluid region)  #TODO ??? how does this work for distance set in MainKratos.py
            for node in self.main_model_part.Nodes:
                distance_value = node.GetSolutionStepValue(KratosMultiphysics.DISTANCE)
                node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, -distance_value)

    def __SetUpInterfaceUtility(self):
        # Create the boundary elements and MLS basis
        settings = KratosMultiphysics.Parameters("""{}""")
        settings.AddEmptyValue("model_part_name").SetString(self.main_model_part.Name + "." + self.GetComputingModelPart().Name)
        settings.AddEmptyValue("boundary_sub_model_part_name").SetString("shifted_boundary")
        settings.AddEmptyValue("conforming_basis").SetBool(self.settings["formulation"]["conforming_basis"].GetBool())
        settings.AddEmptyValue("extension_operator_type").SetString(self.settings["formulation"]["extension_operator_type"].GetString())
        settings.AddEmptyValue("mls_extension_operator_order").SetInt(self.settings["formulation"]["mls_extension_operator_order"].GetInt())
        settings.AddEmptyValue("sbm_interface_condition_name").SetString(self.sbm_interface_condition_name)

        if self.level_set_type == "point-based":
            # Calculate the required neighbors
            elemental_neighbors_process = KratosMultiphysics.GenericFindElementalNeighboursProcess(self.main_model_part)
            elemental_neighbors_process.Execute()

            # Add more specific settings for point-based sbm utility
            settings.AddEmptyValue("skin_model_part_name").SetString("Skin")

            # Store names and interface utilities for all skin model parts
            skin_model_part_names = ["Skin"]
            #skin_model_part_names = ["Cylinder", "VerticalPlate1", "VerticalPlate2", "VerticalPlate3"]  # ["Cylinder", "VerticalPlate1", "VerticalPlate2", "VerticalPlate3"]
            sbm_interface_utilities = []

            # Create an interface utility for all skin model parts
            for skin_model_part_name in skin_model_part_names:
                # Adapt settings
                settings["skin_model_part_name"].SetString(skin_model_part_name)
                if skin_model_part_name == "Cylinder":
                    #TODO NEXT new "find_enclosed_volumes" process instead
                    # to determine enclosed volumes bound by inactive elements and set p=0 for one node if no node's pressure is fixed yet
                    settings.AddEmptyValue("enclosed_area").SetString("negative")
                    #settings.AddEmptyValue("enclosed_area").SetString("none")
                else:
                    settings.AddEmptyValue("enclosed_area").SetString("none")

                # Create interface utility
                sbm_interface_utility = KratosMultiphysics.ShiftedBoundaryPointBasedInterfaceUtility(self.model, settings)
                sbm_interface_utilities.append(sbm_interface_utility)
                KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "New shifted-boundary point-based interface utility created.")

            # Interface flags should be reset for the volume/ computing model part once before skin model parts are (newly) embedded
            sbm_interface_utilities[0].ResetFlags()

            # Locate skin model part points in the volume model part elements and set flags for all skin model parts
            for sbm_interface_utility in sbm_interface_utilities:
                sbm_interface_utility.SetTessellatedBoundaryFlags()
                # To be done after setting tessellated boundary because nodes might be relocated
                sbm_interface_utility.LocateSkinPoints()
                # To be done after locating the skin points because elements in which skin points are located
                # might not be intersected by tessellated skin and might be marked as boundary here
                #sbm_interface_utility.SetInterfaceFlags()

            # Deactivate BOUNDARY elements and nodes which are surrounded by deactivated elements
            #sbm_interface_utilities[0].DeactivateElementsAndNodes()

            # Add Kratos conditions for points at the boundary based on extension operators
            # NOTE that same boundary sub model part is being used here for all skin model parts and their utilities to add conditions
            for i_skin, sbm_interface_utility in enumerate(sbm_interface_utilities):
                sbm_interface_utility.AddSkinIntegrationPointConditions()
                KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Integration point conditions added for skin model part " + skin_model_part_names[i_skin] + ".")

            sbm_interface_utility.SetInterfaceFlags()
            sbm_interface_utilities[0].DeactivateElementsAndNodes()

            # Search for enclosed volumes and fix the pressure of one node if it has not been fixed yet
            #sbm_interface_utilities[0].FixEnclosedVolumesPressure()

        elif self.level_set_type == "discontinuous":
            # Calculate the required neighbors
            elemental_neighbors_process = KratosMultiphysics.GenericFindElementalNeighboursProcess(self.main_model_part)
            elemental_neighbors_process.Execute()

            #
            settings.AddEmptyValue("levelset_variable_name").SetString("ELEMENTAL_DISTANCES")
            sbm_interface_utility = KratosMultiphysics.ShiftedBoundaryMeshlessDiscontinuousInterfaceUtility(self.model, settings)
            KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Shifted-boundary interface utility for discontinuous level-set created.")

            # Add Kratos conditions for points at the boundary based on extension operators
            sbm_interface_utility.CalculateExtensionOperator()

        else:
            # Calculate the required neighbors
            nodal_neighbours_process = KratosMultiphysics.FindGlobalNodalNeighboursProcess(self.main_model_part)
            nodal_neighbours_process.Execute()
            elemental_neighbors_process = KratosMultiphysics.GenericFindElementalNeighboursProcess(self.main_model_part)
            elemental_neighbors_process.Execute()

            settings.AddEmptyValue("levelset_variable_name").SetString("DISTANCE")
            sbm_interface_utility = KratosMultiphysics.ShiftedBoundaryMeshlessInterfaceUtility(self.model, settings)
            KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Shifted-boundary interface utility for continuous level-set created.")

            # Add Kratos conditions for points at the boundary based on extension operators
            sbm_interface_utility.CalculateExtensionOperator()

        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Extension operators were calculated and interface conditions added.")

    def GetDistanceModificationProcess(self):
        if not hasattr(self, '_distance_modification_process'):
            self._distance_modification_process = self.__CreateDistanceModificationProcess()
        return self._distance_modification_process

    def __CreateDistanceModificationProcess(self):
        # Set the distance modification settings according to the level set type
        # Note that the distance modification process is applied to the volume model part
        distance_modification_settings = self.settings["distance_modification_settings"]
        distance_modification_settings.ValidateAndAssignDefaults(self.__GetDistanceModificationDefaultSettings(self.level_set_type))
        aux_full_volume_part_name = self.settings["model_part_name"].GetString() + "." + self.settings["volume_model_part_name"].GetString()
        distance_modification_settings["model_part_name"].SetString(aux_full_volume_part_name)
        return KratosCFD.DistanceModificationProcess(self.model, distance_modification_settings)
