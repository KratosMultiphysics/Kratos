# Importing the Kratos Library
import KratosMultiphysics as KM

# Import applications
import KratosMultiphysics.FluidDynamicsApplication as KM_CFD

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
        default_settings = KM.Parameters(r"""{
            "element_type"                 : "shifted_boundary_weakly_compressible_navier_stokes",
            "boundary_model_parts"         : [],
            "enclosed_areas"               : [],
            "slip_length"                  : 1.0e8,
            "penalty_coefficient"          : 10.0,
            "dynamic_tau"                  : 1.0,
            "level_set_type"               : "point-based",
            "extension_operator_type"      : "MLS",
            "mls_extension_operator_order" : 1,
            "postprocess_drag"             : false,
            "postprocess_velocity"         : false,
            "postprocess_pressure"         : false
        }""")
        formulation_settings.ValidateAndAssignDefaults(default_settings)

        self.element_name = "ShiftedBoundaryWeaklyCompressibleNavierStokes"
        self.condition_name = "NavierStokesWallCondition"
        self.sbm_interface_condition_name = "ShiftedBoundaryWallCondition"
        self.boundary_sub_model_part_name = "ShiftedBoundaryConditions"
        self.level_set_type = formulation_settings["level_set_type"].GetString()
        if self.level_set_type != "point-based":
            err_msg = 'Provided level set type is unknown. Available type for MLS-based SBM is \'point-based\'.'
            raise Exception(err_msg)
        self.element_integrates_in_time = True
        self.element_has_nodal_properties = True
        self.historical_nodal_properties_variables_list = [KM.DENSITY]
        self.non_historical_nodal_properties_variables_list = [KM.SOUND_VELOCITY]

        self.process_info_data[KM.DYNAMIC_TAU] = formulation_settings["dynamic_tau"].GetDouble()
        self.process_info_data[KM.PENALTY_COEFFICIENT] = formulation_settings["penalty_coefficient"].GetDouble()
        self.process_info_data[KM_CFD.SLIP_LENGTH] = formulation_settings["slip_length"].GetDouble()

        # Get names of boundary model parts, whether they enclose and area and set post-processing computations for them
        self.skin_model_part_names = formulation_settings["boundary_model_parts"].GetStringArray()
        self.enclosed_areas = formulation_settings["enclosed_areas"].GetStringArray()
        if len(self.enclosed_areas) == 0:
            self.enclosed_areas = ['none'] * len(self.skin_model_part_names)
        else:
            if len(self.enclosed_areas) != len(self.skin_model_part_names):
                err_msg = 'Provided array of \'enclosed_areas\' must have the same length as \'boundary_model_parts\'.'
                raise Exception(err_msg)
            else:
                for area in self.enclosed_areas:
                    print(area)
                    if area not in ["none", "negative", "positive"]:
                        err_msg = 'Provided enclosed area designation is unknown. Available designations are \'none\', \'negative\' and \'positive\'.'
                        raise Exception(err_msg)

        self.postprocess_drag = formulation_settings["postprocess_drag"].GetBool()
        self.postprocess_velocity = formulation_settings["postprocess_velocity"].GetBool()
        self.postprocess_pressure = formulation_settings["postprocess_pressure"].GetBool()


def CreateSolver(model, custom_settings):
    return NavierStokesShiftedBoundaryMonolithicSolver(model, custom_settings)

class NavierStokesShiftedBoundaryMonolithicSolver(FluidSolver):

    @classmethod
    def GetDefaultParameters(cls):
        ##settings string in json format
        default_settings = KM.Parameters("""
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
            "maximum_iterations": 7,
            "echo_level": 0,
            "time_order": 2,
            "time_scheme": "bdf2",
            "compute_reactions": false,
            "analysis_type": "non_linear",
            "reform_dofs_at_each_step": false,
            "consider_periodic_conditions": false,
            "assign_neighbour_elements_to_conditions": true,
            "relative_velocity_tolerance": 1e-5,
            "absolute_velocity_tolerance": 1e-7,
            "relative_pressure_tolerance": 1e-5,
            "absolute_pressure_tolerance": 1e-7,
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
        self.skin_model_part_names = self.shifted_boundary_formulation.skin_model_part_names
        self.enclosed_areas = self.shifted_boundary_formulation.enclosed_areas
        self.postprocess_drag = self.shifted_boundary_formulation.postprocess_drag
        self.postprocess_velocity = self.shifted_boundary_formulation.postprocess_velocity
        self.postprocess_pressure = self.shifted_boundary_formulation.postprocess_pressure
        self.boundary_sub_model_part_name = self.shifted_boundary_formulation.boundary_sub_model_part_name

        # Create a skin model part
        for skin_mp_name in self.skin_model_part_names:
            self.model.CreateModelPart(skin_mp_name)

        # Create a model part for skin integration points
        self.skin_point_model_part = self.model.CreateModelPart("SkinPoints")

        KM.Logger.PrintInfo(self.__class__.__name__, "Construction of NavierStokesShiftedBoundaryMonolithicSolver finished.")

    def AddVariables(self):
        self.main_model_part.AddNodalSolutionStepVariable(KM.PRESSURE)
        self.main_model_part.AddNodalSolutionStepVariable(KM.VELOCITY)
        self.main_model_part.AddNodalSolutionStepVariable(KM.ACCELERATION)
        self.main_model_part.AddNodalSolutionStepVariable(KM.MESH_VELOCITY)
        self.main_model_part.AddNodalSolutionStepVariable(KM.BODY_FORCE)
        self.main_model_part.AddNodalSolutionStepVariable(KM.NODAL_AREA)
        self.main_model_part.AddNodalSolutionStepVariable(KM.REACTION)
        self.main_model_part.AddNodalSolutionStepVariable(KM.REACTION_WATER_PRESSURE)
        self.main_model_part.AddNodalSolutionStepVariable(KM.NORMAL)
        self.main_model_part.AddNodalSolutionStepVariable(KM.EXTERNAL_PRESSURE)
        self.main_model_part.AddNodalSolutionStepVariable(KM.DISTANCE)              # Use distance variable for node relocation and voting on positive/ negative side
        self.main_model_part.AddNodalSolutionStepVariable(KM_CFD.EMBEDDED_WET_PRESSURE)          # Post-process variable (stores the fluid nodes pressure and is set to 0 in the structure ones)
        self.main_model_part.AddNodalSolutionStepVariable(KM_CFD.EMBEDDED_WET_VELOCITY)          # Post-process variable (stores the fluid nodes velocity and is set to 0 in the structure ones)

        # Adding variables required for the nodal material properties
        if self.element_has_nodal_properties:
            for variable in self.historical_nodal_properties_variables_list:
                self.main_model_part.AddNodalSolutionStepVariable(variable)

        # Add nodal variables to skin model part
        for skin_mp_name in self.skin_model_part_names:
            skin_mp = self.model.GetModelPart(skin_mp_name)
            skin_mp.AddNodalSolutionStepVariable(KM.VELOCITY)
            skin_mp.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
            skin_mp.AddNodalSolutionStepVariable(KM.POSITIVE_FACE_PRESSURE)
            skin_mp.AddNodalSolutionStepVariable(KM.NEGATIVE_FACE_PRESSURE)
            skin_mp.AddNodalSolutionStepVariable(KM.POSITIVE_FACE_FLUID_VELOCITY)
            skin_mp.AddNodalSolutionStepVariable(KM.NEGATIVE_FACE_FLUID_VELOCITY)

        # Add nodal variables for skin integration points model part
        self.skin_point_model_part.AddNodalSolutionStepVariable(KM.POSITIVE_FACE_PRESSURE)
        self.skin_point_model_part.AddNodalSolutionStepVariable(KM.NEGATIVE_FACE_PRESSURE)
        self.skin_point_model_part.AddNodalSolutionStepVariable(KM.POSITIVE_FACE_FLUID_VELOCITY)
        self.skin_point_model_part.AddNodalSolutionStepVariable(KM.NEGATIVE_FACE_FLUID_VELOCITY)
        self.skin_point_model_part.AddNodalSolutionStepVariable(KM.TRACTION_FROM_FLUID_PRESSURE)
        self.skin_point_model_part.AddNodalSolutionStepVariable(KM.TRACTION_FROM_FLUID_STRESS)

        KM.Logger.PrintInfo(self.__class__.__name__, "Shifted-boundary fluid solver variables added correctly.")

    def ImportModelPart(self):
        super(NavierStokesShiftedBoundaryMonolithicSolver, self).ImportModelPart()

    def PrepareModelPart(self):
        if self.main_model_part.ProcessInfo[KM.IS_RESTARTED]:
            # Delete conditions related to interface utility - the interface utility and conditions will be recreated
            boundary_sub_model_part = self.GetComputingModelPart().GetSubModelPart(self.boundary_sub_model_part_name)
            KM.VariableUtils().SetFlag(KM.TO_ERASE, True, boundary_sub_model_part.Conditions)
            boundary_sub_model_part.RemoveConditions(KM.TO_ERASE)

        # Call the fluid solver PrepareModelPart()
        super(NavierStokesShiftedBoundaryMonolithicSolver, self).PrepareModelPart()

        # Set the extra requirements of the shifted-boundary formulation
        if not self.main_model_part.ProcessInfo[KM.IS_RESTARTED]:
            # Set the shifted-boundary formulation configuration
            self.__SetShiftedBoundaryFormulation()

        # Clone the solution step data for skin and skin points model parts
        t =  self.GetComputingModelPart().ProcessInfo[KM.TIME]
        step = self.GetComputingModelPart().ProcessInfo[KM.STEP]
        for skin_mp_name in self.skin_model_part_names:
            skin_mp = self.model.GetModelPart(skin_mp_name)
            skin_mp.CloneTimeStep(t)
            skin_mp.ProcessInfo[KM.STEP] = step
        self.skin_point_model_part.CloneTimeStep(t)
        self.skin_point_model_part.ProcessInfo[KM.STEP] = step

    def Initialize(self):
        # If the solver requires an instance of the stabilized shifted boundary formulation class, set the process info variables
        if hasattr(self, 'shifted_boundary_formulation'):
            self.shifted_boundary_formulation.SetProcessInfo(self.GetComputingModelPart())

        # Construct and initialize the solution strategy
        # NOTE "Error: Constitutive Law not initialized for Element ShiftedBoundaryFluidElement #105033"
        # if strategy is initialized after set up of interface utility (deactivation of elements?)
        solution_strategy = self._GetSolutionStrategy()
        solution_strategy.SetEchoLevel(self.settings["echo_level"].GetInt())
        solution_strategy.Initialize()

        # Deactivate darcy term  #TODO necessary?
        # for ele in self.GetComputingModelPart().Elements:
        #     ele.SetValue(KM_CFD.RESISTANCE, 0.0)

        # Create shifted-boundary meshless interface utility and calculate extension operator requiring nodal and elemental neighbors
        self.__SetUpInterfaceUtility()

        KM.Logger.PrintInfo(self.__class__.__name__, "Solver initialization finished.")

    def AdvanceInTime(self, current_time):
        new_time = super(NavierStokesShiftedBoundaryMonolithicSolver, self).AdvanceInTime(current_time)

        # Clone the solution step data for skin and skin points model parts
        step = self.GetComputingModelPart().ProcessInfo[KM.STEP]
        for skin_mp_name in self.skin_model_part_names:
            skin_mp = self.model.GetModelPart(skin_mp_name)
            skin_mp.CloneTimeStep(new_time)
            skin_mp.ProcessInfo[KM.STEP] = step
        self.skin_point_model_part.CloneTimeStep(new_time)
        self.skin_point_model_part.ProcessInfo[KM.STEP] = step

        return new_time

    def InitializeSolutionStep(self):
        # Compute the BDF coefficients
        (self.time_discretization).ComputeAndSaveBDFCoefficients(self.GetComputingModelPart().ProcessInfo)

        # Call the base solver InitializeSolutionStep()
        super(NavierStokesShiftedBoundaryMonolithicSolver, self).InitializeSolutionStep()

    def FinalizeSolutionStep(self):
        # Compute Variables for skin model parts in sbm utilities
        for sbm_interface_utility in self.sbm_interface_utilities:
            if self.postprocess_drag:
                sbm_interface_utility.CalculateSkinDrag()
            if self.postprocess_velocity:
                sbm_interface_utility.CalculateVelocityAtSkinNodes()
            if self.postprocess_pressure:
                sbm_interface_utility.CalculatePressureAtSkinNodes()

        # Call the base solver FinalizeSolutionStep()
        super(NavierStokesShiftedBoundaryMonolithicSolver, self).FinalizeSolutionStep()

    def _SetNodalProperties(self):
        set_density = KM.DENSITY in self.historical_nodal_properties_variables_list
        set_sound_velocity = KM.SOUND_VELOCITY in self.non_historical_nodal_properties_variables_list

        # Get density and dynamic viscosity from the properties of the first element
        for ele in self.main_model_part.Elements:
            # Get DENSITY from properties
            if set_density:
                rho = ele.Properties.GetValue(KM.DENSITY)
                if rho <= 0.0:
                    raise Exception("DENSITY set to {0} in Properties {1}, positive number expected.".format(rho,ele.Properties.Id))
            # Get SOUND_VELOCITY
            if set_sound_velocity:
                if ele.Properties.Has(KM.SOUND_VELOCITY):
                    sound_velocity = ele.Properties.GetValue(KM.SOUND_VELOCITY)
                else:
                    # NOTE this is the default sound velocity value to deactivate compressibility term of elements
                    sound_velocity = 1.0e+12
                    KM.Logger.PrintWarning('No \'SOUND_VELOCITY\' value found in Properties {0}. Setting default value {1}'.format(ele.Properties.Id, sound_velocity))
                if sound_velocity <= 0.0:
                    raise Exception("SOUND_VELOCITY set to {0} in Properties {1}, positive number expected.".format(sound_velocity, ele.Properties.Id))
            break
        else:
            raise Exception("No fluid elements found in the main model part.")

        # Transfer the obtained properties to the nodes
        if set_density:
            KM.VariableUtils().SetVariable(KM.DENSITY, rho, self.main_model_part.Nodes)
        if set_sound_velocity:
            KM.VariableUtils().SetNonHistoricalVariable(KM.SOUND_VELOCITY, sound_velocity, self.main_model_part.Nodes)

    def __SetShiftedBoundaryFormulation(self):
        # Save the formulation settings in the ProcessInfo
        self.shifted_boundary_formulation.SetProcessInfo(self.main_model_part)

    def __SetUpInterfaceUtility(self):
        # Create the boundary elements and MLS basis
        settings = KM.Parameters("""{}""")
        settings.AddEmptyValue("model_part_name").SetString(self.main_model_part.Name + "." + self.GetComputingModelPart().Name)
        settings.AddEmptyValue("boundary_sub_model_part_name").SetString(self.boundary_sub_model_part_name)
        settings.AddEmptyValue("extension_operator_type").SetString(self.settings["formulation"]["extension_operator_type"].GetString())
        settings.AddEmptyValue("mls_extension_operator_order").SetInt(self.settings["formulation"]["mls_extension_operator_order"].GetInt())
        settings.AddEmptyValue("sbm_interface_condition_name").SetString(self.sbm_interface_condition_name)

        #if self.level_set_type == "point-based":
        #n_dim = self.main_model_part.ProcessInfo[KM.DOMAIN_SIZE]

        # Calculate the required neighbors
        elemental_neighbors_process = KM.GenericFindElementalNeighboursProcess(self.main_model_part)
        elemental_neighbors_process.Execute()

        # Create an interface utility for all skin model part names
        self.sbm_interface_utilities = []
        for skin_model_part_name, enclosed_area in zip(self.skin_model_part_names, self.enclosed_areas):
            # Adapt settings
            settings.AddEmptyValue("skin_model_part_name").SetString(skin_model_part_name)
            settings.AddEmptyValue("enclosed_area").SetString(enclosed_area)

            # Create interface utility
            sbm_interface_utility = KM.ShiftedBoundaryPointBasedInterfaceUtility(self.model, settings)
            self.sbm_interface_utilities.append(sbm_interface_utility)
            KM.Logger.PrintInfo(self.__class__.__name__, "New shifted-boundary point-based interface utility created for skin model part '" + skin_model_part_name + "'.")

        if len(self.sbm_interface_utilities) == 1:
            self.sbm_interface_utilities[0].CalculateAndAddPointBasedInterface()
            KM.Logger.PrintInfo(self.__class__.__name__, "Extension operators were calculated and interface conditions added.")

        elif len(self.sbm_interface_utilities) > 1:
            # Interface flags should be reset for the volume/ computing model part once before skin model parts are (newly) embedded
            self.sbm_interface_utilities[0].ResetFlags()

            # Set boundary flags and locate skin model part points in the volume model part elements for all skin model parts
            for sbm_interface_utility in self.sbm_interface_utilities:
                sbm_interface_utility.SetTessellatedBoundaryFlagsAndRelocateSmallDistanceNodes()
                # To be done after setting tessellated boundary because nodes might be relocated
                sbm_interface_utility.LocateSkinPoints()

            # To be done after locating the skin points because elements in which skin points are located
            # might not be intersected by tessellated skin and might be marked as boundary here
            self.sbm_interface_utilities[0].SetInterfaceFlags()

            # Deactivate BOUNDARY elements and nodes which are surrounded by deactivated elements
            self.sbm_interface_utilities[0].DeactivateElementsAndNodes()

            # Add Kratos conditions for points at the boundary based on extension operators
            # NOTE that the same boundary sub model part is being used here for all skin model parts and their utilities to add conditions
            for i_skin, sbm_interface_utility in enumerate(self.sbm_interface_utilities):
                sbm_interface_utility.CalculateAndAddSkinIntegrationPointConditions()
                KM.Logger.PrintInfo(self.__class__.__name__, "Integration point conditions added for skin model part '" + self.skin_model_part_names[i_skin] + "'.")

            #TODO Search for enclosed volumes and fix the pressure of one node if it has not been fixed yet? (instead of defining enclosed_areas)
            #sbm_interface_utilities[0].FixEnclosedVolumesPressure()

            KM.Logger.PrintInfo(self.__class__.__name__, "Extension operators were calculated and interface conditions added.")

        else:
            KM.Logger.PrintWarning('Shifted-boundary interface utility was not set up because no boundary model part was given.')

