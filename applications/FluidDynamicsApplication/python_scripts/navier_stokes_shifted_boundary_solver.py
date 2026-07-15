# Importing the Kratos Library
import KratosMultiphysics as KM

# Import Kratos utilities
from KratosMultiphysics.kratos_utilities import CheckIfApplicationsAvailable

# Import applications
import KratosMultiphysics.FluidDynamicsApplication as KM_CFD
mesh_moving_available = CheckIfApplicationsAvailable("MeshMovingApplication")
if mesh_moving_available:
    import KratosMultiphysics.MeshMovingApplication as KM_MeshMoving

# Import base class file
from KratosMultiphysics.FluidDynamicsApplication.fluid_solver import FluidSolver
from KratosMultiphysics.FluidDynamicsApplication import check_and_prepare_model_process_fluid

import datetime
import numpy

#from bin.RelWithDebInfo import KratosMultiphysics


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
            "element_type"                   : "shifted_boundary_weakly_compressible_navier_stokes",
            "boundary_model_parts"           : [],
            "enclosed_areas"                 : [],
            "deactivate_unstable_clusters"   : false,
            "slip_length"                    : 1.0e8,
            "charact_length"                 : 1.0,
            "penalty_coefficient"            : 10.0,
            "penalty_coefficient_tangential" : 10.0,
            "dynamic_tau"                    : 1.0,
            "level_set_type"                 : "point-based",
            "extension_operator_type"        : "MLS",
            "mls_extension_operator_order"   : 1,
            "postprocess_skin_points"        : false,
            "postprocess_skin_nodes"         : false
        }""")
        formulation_settings.ValidateAndAssignDefaults(default_settings)

        self.element_name = "ShiftedBoundaryWeaklyCompressibleNavierStokes"
        self.condition_name = "NavierStokesWallCondition"
        self.boundary_wall_condition_name = "ShiftedBoundaryWallCondition"
        self.boundary_sub_model_part_name = "ShiftedBoundaryConditions"
        self.level_set_type = formulation_settings["level_set_type"].GetString()
        if self.level_set_type != "point-based":
            err_msg = 'Provided level set type is unknown. Only available type for Navier-Stokes SBM is \'point-based\'.'
            raise Exception(err_msg)
        self.element_integrates_in_time = True
        self.element_has_nodal_properties = True
        self.historical_nodal_properties_variables_list = [KM.DENSITY]
        self.non_historical_nodal_properties_variables_list = [KM.SOUND_VELOCITY]

        self.process_info_data[KM.DYNAMIC_TAU] = formulation_settings["dynamic_tau"].GetDouble()
        self.process_info_data[KM.PENALTY_COEFFICIENT] = formulation_settings["penalty_coefficient"].GetDouble()
        self.process_info_data[KM.PENALTY_COEFFICIENT_TANGENTIAL] = formulation_settings["penalty_coefficient_tangential"].GetDouble()
        self.process_info_data[KM_CFD.SLIP_LENGTH] = formulation_settings["slip_length"].GetDouble()
        self.process_info_data[KM_CFD.EMBEDDED_CHARACT_LENGTH] = formulation_settings["charact_length"].GetDouble()

        # Get names of boundary model parts, whether they enclose an area
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
                    if area not in ["none", "negative", "positive"]:
                        err_msg = 'Provided enclosed area designation is unknown. Available designations are \'none\', \'negative\' and \'positive\'.'
                        raise Exception(err_msg)

        # Decide whether enclosed volumes of the fluid domain created by embedded geometries should be deactivated if no degree of freedom is fixed (unstable)
        self.deactivate_unstable_clusters = formulation_settings["deactivate_unstable_clusters"].GetBool()

        # (De-)activate post-processing routines for skin integration points and nodes
        self.postprocess_skin_points = formulation_settings["postprocess_skin_points"].GetBool()
        self.postprocess_skin_nodes = formulation_settings["postprocess_skin_nodes"].GetBool()


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
            "enforce_element_and_conditions_replacement": true,
            "material_import_settings": {
                "materials_filename": ""
            },
            "maximum_iterations": 7,
            "echo_level": 0,
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
            },
            "fm_ale_settings": {
                "fm_ale_step_frequency": 0,
                "mesh_movement": "implicit",
                "embedded_velocity_calculation": "from_fluid_mesh_velocity",
                "rbf_interpolation_search_radius": 0.0,
                "rbf_interpolation_edge_factor": 2.0,
                "fm_ale_solver_settings": {
                }
            }
        }""")

        default_settings.AddMissingParameters(super(NavierStokesShiftedBoundaryMonolithicSolver, cls).GetDefaultParameters())
        return default_settings

    @classmethod
    def _GetFmAleImplicitDefaultSettings(cls):
        return KM.Parameters("""
        {
            "virtual_model_part_name": "VirtualModelPart",
            "structure_model_part_name": "",
            "linear_solver_settings": {
                "solver_type": "cg",
                "tolerance": 1.0e-8,
                "max_iteration": 1000
            },
            "embedded_nodal_variable_settings": {
                "gradient_penalty_coefficient": 1.0e-3,
                "linear_solver_settings": {
                    "preconditioner_type": "amg",
                    "solver_type": "amgcl",
                    "smoother_type": "ilu0",
                    "krylov_type": "cg",
                    "max_iteration": 1000,
                    "verbosity": 0,
                    "tolerance": 1e-8,
                    "scaling": false,
                    "block_size": 1,
                    "use_block_matrices_if_possible": true
                }
            }
        }
        """)

    def ValidateSettings(self):
        """ Overriding python_solver ValidateSettings to validate the FM-ALE settings
        """
        super(NavierStokesShiftedBoundaryMonolithicSolver, self).ValidateSettings()

        self.settings["fm_ale_settings"].ValidateAndAssignDefaults(self.GetDefaultParameters()["fm_ale_settings"])
        if self.settings["fm_ale_settings"]["fm_ale_step_frequency"].GetInt() > 0:
            mesh_movement = self.settings["fm_ale_settings"]["mesh_movement"].GetString()
            if mesh_movement != "implicit":
                raise Exception("Provided mesh movement \'" + mesh_movement + "\'. Available options are \'implicit\'.")
            self.settings["fm_ale_settings"]["fm_ale_solver_settings"].ValidateAndAssignDefaults(self._GetFmAleImplicitDefaultSettings())

    def __init__(self, model, custom_settings):
        super(NavierStokesShiftedBoundaryMonolithicSolver, self).__init__(model, custom_settings)

        self.min_buffer_size = 3
        self.shifted_boundary_formulation = ShiftedBoundaryFormulation(self.settings["formulation"])
        self.element_name = self.shifted_boundary_formulation.element_name
        self.condition_name = self.shifted_boundary_formulation.condition_name
        self.sbm_wall_condition_name = self.shifted_boundary_formulation.boundary_wall_condition_name + str(self.settings["domain_size"].GetInt()) + "D"
        self.level_set_type = self.shifted_boundary_formulation.level_set_type
        self.element_integrates_in_time = self.shifted_boundary_formulation.element_integrates_in_time
        self.element_has_nodal_properties = self.shifted_boundary_formulation.element_has_nodal_properties
        self.historical_nodal_properties_variables_list = self.shifted_boundary_formulation.historical_nodal_properties_variables_list
        self.non_historical_nodal_properties_variables_list = self.shifted_boundary_formulation.non_historical_nodal_properties_variables_list
        self.skin_model_part_names = self.shifted_boundary_formulation.skin_model_part_names
        self.enclosed_areas = self.shifted_boundary_formulation.enclosed_areas
        self.deactivate_unstable_clusters = self.shifted_boundary_formulation.deactivate_unstable_clusters
        self.postprocess_skin_points = self.shifted_boundary_formulation.postprocess_skin_points
        self.postprocess_skin_nodes = self.shifted_boundary_formulation.postprocess_skin_nodes
        self.boundary_sub_model_part_name = self.shifted_boundary_formulation.boundary_sub_model_part_name

        # Create a skin model part and skin points model part
        for skin_mp_name in self.skin_model_part_names:
            if not self.model.HasModelPart(skin_mp_name):
                skin_mp = self.model.CreateModelPart(skin_mp_name)
            if not self.model.HasModelPart(skin_mp_name + "Points"):
                self.model.CreateModelPart(skin_mp_name + "Points")

        # If the FM-ALE is required, do a first call to __GetFmAleVirtualModelPart
        # NOTE that this will create the virtual model part in the model
        if self._FmAleIsActive():
            self.__GetFmAleVirtualModelPart()

        KM.Logger.PrintInfo(self.__class__.__name__, "Construction of NavierStokesShiftedBoundaryMonolithicSolver finished.")

    def AddVariables(self):
        self.main_model_part.AddNodalSolutionStepVariable(KM.PRESSURE)
        self.main_model_part.AddNodalSolutionStepVariable(KM.VELOCITY)
        self.main_model_part.AddNodalSolutionStepVariable(KM.ACCELERATION)
        self.main_model_part.AddNodalSolutionStepVariable(KM.MESH_VELOCITY)
        self.main_model_part.AddNodalSolutionStepVariable(KM.BODY_FORCE)  # used?
        self.main_model_part.AddNodalSolutionStepVariable(KM.NODAL_AREA)  # used?
        self.main_model_part.AddNodalSolutionStepVariable(KM.REACTION)  # used?
        self.main_model_part.AddNodalSolutionStepVariable(KM.REACTION_WATER_PRESSURE)  # used?
        self.main_model_part.AddNodalSolutionStepVariable(KM.NORMAL)  # used?
        self.main_model_part.AddNodalSolutionStepVariable(KM.EXTERNAL_PRESSURE)  # used?
        self.main_model_part.AddNodalSolutionStepVariable(KM.DISTANCE)                       # Use distance variable for voting on positive/ negative side
        #self.main_model_part.AddNodalSolutionStepVariable(KM_CFD.EMBEDDED_WET_PRESSURE)     # Post-process variable (stores the fluid nodes pressure and is set to 0 in the structure ones)
        #self.main_model_part.AddNodalSolutionStepVariable(KM_CFD.EMBEDDED_WET_VELOCITY)     # Post-process variable (stores the fluid nodes velocity and is set to 0 in the structure ones)
        #self.main_model_part.AddNodalSolutionStepVariable(KM.EMBEDDED_VELOCITY)

        # Add variables required for the FM-ALE algorithm
        if self._FmAleIsActive():
            self.main_model_part.AddNodalSolutionStepVariable(KM.MESH_DISPLACEMENT)
            self.main_model_part.AddNodalSolutionStepVariable(KM.MESH_REACTION)

        # Add variables required for the nodal material properties
        if self.element_has_nodal_properties:
            for variable in self.historical_nodal_properties_variables_list:
                self.main_model_part.AddNodalSolutionStepVariable(variable)

        # Add nodal variables to skin model part and skin integration points model part
        for skin_mp_name in self.skin_model_part_names:
            skin_mp = self.model.GetModelPart(skin_mp_name)
            skin_mp.AddNodalSolutionStepVariable(KM.VELOCITY)
            skin_mp.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
            skin_mp.AddNodalSolutionStepVariable(KM.POSITIVE_FACE_PRESSURE)
            skin_mp.AddNodalSolutionStepVariable(KM.NEGATIVE_FACE_PRESSURE)
            skin_mp.AddNodalSolutionStepVariable(KM.POSITIVE_FACE_FLUID_VELOCITY)
            skin_mp.AddNodalSolutionStepVariable(KM.NEGATIVE_FACE_FLUID_VELOCITY)
            skin_mp_points = self.model.GetModelPart(skin_mp_name + "Points")
            skin_mp_points.AddNodalSolutionStepVariable(KM.NORMAL)
            skin_mp_points.AddNodalSolutionStepVariable(KM.POSITIVE_FACE_PRESSURE)
            skin_mp_points.AddNodalSolutionStepVariable(KM.NEGATIVE_FACE_PRESSURE)
            skin_mp_points.AddNodalSolutionStepVariable(KM.POSITIVE_FACE_FLUID_VELOCITY)
            skin_mp_points.AddNodalSolutionStepVariable(KM.NEGATIVE_FACE_FLUID_VELOCITY)
            skin_mp_points.AddNodalSolutionStepVariable(KM.TRACTION_FROM_FLUID_PRESSURE)
            skin_mp_points.AddNodalSolutionStepVariable(KM.TRACTION_FROM_FLUID_STRESS)
            skin_mp_points.AddNodalSolutionStepVariable(KM.DRAG_FORCE)

        KM.Logger.PrintInfo(self.__class__.__name__, "Shifted-boundary fluid solver variables added correctly.")

    def AddDofs(self):
        # Add formulation DOFs and reactions
        super(NavierStokesShiftedBoundaryMonolithicSolver, self).AddDofs()

        # Add mesh motion problem DOFs for the FM-ALE algorithm
        if self._FmAleIsActive():
            dofs_and_reactions_to_add = []
            dofs_and_reactions_to_add.append(["MESH_DISPLACEMENT_X", "MESH_REACTION_X"])
            dofs_and_reactions_to_add.append(["MESH_DISPLACEMENT_Y", "MESH_REACTION_Y"])
            dofs_and_reactions_to_add.append(["MESH_DISPLACEMENT_Z", "MESH_REACTION_Z"])
            KM.VariableUtils.AddDofsList(dofs_and_reactions_to_add, self.main_model_part)

            KM.Logger.PrintInfo(self.__class__.__name__, "FM-ALE DOFs added correctly.")

    def PrepareModelPart(self):
        # Delete SBM conditions if the simulation is restarted - the shifted-boundary utility and conditions will be recreated
        if self.main_model_part.ProcessInfo[KM.IS_RESTARTED]:
            boundary_sub_model_part = self.GetComputingModelPart().GetSubModelPart(self.boundary_sub_model_part_name)
            for cond in boundary_sub_model_part.Conditions:
                self.main_model_part.GetCondition(cond.Id).Set(KM.TO_ERASE, True)
            self.main_model_part.RemoveConditions(KM.TO_ERASE)

        # Call the fluid solver PrepareModelPart()
        #NOTE this creates the computational fluid model part (if not restarted)
        super(NavierStokesShiftedBoundaryMonolithicSolver, self).PrepareModelPart()

        # Set the extra requirements of the shifted-boundary formulation
        if not self.main_model_part.ProcessInfo[KM.IS_RESTARTED]:
            # Create sub model part for SBM conditions
            if not self.GetComputingModelPart().HasSubModelPart(self.boundary_sub_model_part_name):
                self.GetComputingModelPart().CreateSubModelPart(self.boundary_sub_model_part_name)

            # Set the shifted-boundary formulation configuration
            self.__SetShiftedBoundaryFormulation()

        # Create shifted-boundary utility
        self.__CreateShiftedBoundaryUtilities()
        # Flag BOUNDARY elements for calculating the metric for an initial remeshing (MMG)
        #self.__FlagBoundaryElements()

        # Clone the solution step data for skin and skin points model parts
        t =  self.GetComputingModelPart().ProcessInfo[KM.TIME]
        step = self.GetComputingModelPart().ProcessInfo[KM.STEP]
        for skin_mp_name in self.skin_model_part_names:
            skin_mp = self.model.GetModelPart(skin_mp_name)
            skin_mp.CloneTimeStep(t)
            skin_mp.ProcessInfo[KM.STEP] = step
            skin_mp_points = self.model.GetModelPart(skin_mp_name+"Points")
            skin_mp_points.CloneTimeStep(t)
            skin_mp_points.ProcessInfo[KM.STEP] = step

    def Initialize(self):
        # Run check and prepare process again
        #NOTE that this is necessary after 'initial_remeshing' of a MMG process to get correct parent elements for wall boundary conditions
        # prepare_model_part_settings = KM.Parameters("{}")
        # prepare_model_part_settings.AddValue("volume_model_part_name",self.settings["volume_model_part_name"])
        # prepare_model_part_settings.AddValue("skin_parts",self.settings["skin_parts"])
        # prepare_model_part_settings.AddValue("assign_neighbour_elements_to_conditions",self.settings["assign_neighbour_elements_to_conditions"])
        # check_and_prepare_model_process_fluid.CheckAndPrepareModelProcessFluid(self.main_model_part, prepare_model_part_settings).Execute()

        # If the solver requires an instance of the stabilized shifted boundary formulation class, set the process info variables
        #TODO not necessary?
        # if hasattr(self, 'shifted_boundary_formulation'):
        #    self.shifted_boundary_formulation.SetProcessInfo(self.GetComputingModelPart())

        # Construct and initialize the solution strategy
        # NOTE This needs to be done before immersing the shifted boundary, otherwise there is a constitutive law error b/c of deactivation of elements(?).
        solution_strategy = self._GetSolutionStrategy()
        solution_strategy.SetEchoLevel(self.settings["echo_level"].GetInt())
        solution_strategy.Initialize()

        #TODO would this be useful?
        # Deactivate darcy term
        # for ele in self.GetComputingModelPart().Elements:
        #     ele.SetValue(KM_CFD.RESISTANCE, 0.0)

        # Create skin points in skin points model part from the discretized skin model part if the skin points model part is empty and the skin model part is not empty.
        for skin_mp_name in self.skin_model_part_names:
            skin_mp = self.model.GetModelPart(skin_mp_name)
            skin_mp_points = self.model.GetModelPart(skin_mp_name+"Points")
            if skin_mp_points.NumberOfNodes() == 0 and skin_mp.NumberOfNodes() != 0:
                self.__FillSkinPointsFromDiscModelPart(skin_mp_name, skin_mp, skin_mp_points)

        # If required, initialize the FM-ALE utility
        if self._FmAleIsActive():
            self.fm_ale_step = 1
            # Fill the virtual model part geometry. Note that the mesh moving util is created in this first call
            self.__GetFmAleUtility().Initialize(self.main_model_part)

            #TODO REMOVE
            self.main_model_part.ProcessInfo.SetValue(KM.SOUND_VELOCITY, 1e12)
            self.GetComputingModelPart().ProcessInfo.SetValue(KM.SOUND_VELOCITY, 1e12)

        #else:
            # Call shifted-boundary utility methods to immerse the boundary and calculate extension operators requiring nodal and elemental neighbors
        #TODO why is this necessary??? otherwise SegFault bc of DOF set??
        self.__ImmerseShiftedBoundaries()

        KM.Logger.PrintInfo(self.__class__.__name__, "Solver initialization finished.")

    def AdvanceInTime(self, current_time):
        new_time = super(NavierStokesShiftedBoundaryMonolithicSolver, self).AdvanceInTime(current_time)

        # Clone the solution step data for skin and skin points model parts
        step = self.GetComputingModelPart().ProcessInfo[KM.STEP]
        for skin_mp_name in self.skin_model_part_names:
            skin_mp = self.model.GetModelPart(skin_mp_name)
            skin_mp.CloneTimeStep(new_time)
            skin_mp.ProcessInfo[KM.STEP] = step
            skin_mp_points = self.model.GetModelPart(skin_mp_name+"Points")
            skin_mp_points.CloneTimeStep(new_time)
            skin_mp_points.ProcessInfo[KM.STEP] = step

        # Advance the FM-ALE virtual model part
        if self._FmAleIsActive():
            self.__GetFmAleVirtualModelPart().ProcessInfo[KM.STEP] += 1
            self.__GetFmAleVirtualModelPart().ProcessInfo[KM.TIME] = new_time

        return new_time

    def InitializeSolutionStep(self):
        # Compute the BDF coefficients
        (self.time_discretization).ComputeAndSaveBDFCoefficients(self.GetComputingModelPart().ProcessInfo)

        if self.__IsFmAleStep():
            # Set the virtual mesh values (VELOCITY and PRESSURE) from the background fluid mesh
            self.__GetFmAleUtility().SetVirtualMeshValuesFromOriginMesh()

        # Call the base solver InitializeSolutionStep()
        # NOTE that here the system is resized if demanded
        # TODO deactivate DOFs before and ReformDofSetAtEachStep or keep all DOFs for moving boundary and add BCs??
        #super(NavierStokesShiftedBoundaryMonolithicSolver, self).InitializeSolutionStep()

        if self.__IsFmAleStep():
            # Perform the FM-ALE operations
            self.__DoFmAleMeshMovementAndProjection()

            # Update the skin points model parts from the skin model parts and set the EMBEDDED_VELOCITY in the skin points
            self.__UpdateSkinPointsFromDiscModelParts()

            # Call shifted-boundary utility methods to immerse the updated boundary and calculate extension operators
            self.__ReImmerseShiftedBoundaries()

        super(NavierStokesShiftedBoundaryMonolithicSolver, self).InitializeSolutionStep()

    def SolveSolutionStep(self):
        # TODO boundary movement and shifted-boundaries immersion need to be done here in case of non-linear coupling iterations??
        # problematic for strategy initialization??
        # if self.__IsFmAleStep():
        #     # Perform the FM-ALE operations
        #     self.__DoFmAleMeshMovementAndProjection()

        #     # Calculate the EMBEDDED_VELOCITY
        #     self.__CalculateEmbeddedVelocity()
        #     # TODO Update the skin points model parts from the skin model parts

        #     # Call shifted-boundary utility methods to immerse the updated boundary and calculate extension operators
        #     self.__ReImmerseShiftedBoundaries()

        #     #super(NavierStokesShiftedBoundaryMonolithicSolver, self).InitializeSolutionStep()

        # Call the base SolveSolutionStep to solve the embedded CFD problem
        is_converged = super(NavierStokesShiftedBoundaryMonolithicSolver, self).SolveSolutionStep()

        # # Undo the FM-ALE virtual mesh movement
        # if self._FmAleIsActive():
        #     self.__GetFmAleUtility().UndoMeshMovement()

        return is_converged

    def FinalizeSolutionStep(self):
        # Compute Variables for skin model parts in sbm utilities
        #TODO only call these functions if it is an output step?! AND call CalculateVariablesAtSkinPoints() for every FSI iteration?
        time_prev = datetime.datetime.now().replace(microsecond=0)
        if self.postprocess_skin_nodes:
            for sbm_utility in self.sbm_utilities:
                sbm_utility.CalculateVariablesAtSkinPointsAndNodes()  # This method calculates the variables at the skin points as well
            self.__PrintAndResetTimer(time_prev, "post-process skin variables")
        elif self.postprocess_skin_points:
            for sbm_utility in self.sbm_utilities:
                sbm_utility.CalculateVariablesAtSkinPoints()
            self.__PrintAndResetTimer(time_prev, "post-process skin variables")

        if self._FmAleIsActive():
            # Undo the FM-ALE virtual mesh movement
            self.__GetFmAleUtility().UndoMeshMovement()
            # Reset or update FM-ALE step counter
            if (self.__IsFmAleStep()):
                self.fm_ale_step = 1
            else:
                self.fm_ale_step += 1

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

    def __CreateShiftedBoundaryUtilities(self):
        # Create the boundary elements and MLS parameter basis
        settings = KM.Parameters("""{}""")
        settings.AddEmptyValue("model_part_name").SetString(self.main_model_part.Name)  #TODO + "." + self.GetComputingModelPart().Name)
        settings.AddEmptyValue("boundary_sub_model_part_name").SetString(self.main_model_part.Name + "." + self.GetComputingModelPart().Name + "." + self.boundary_sub_model_part_name)
        settings.AddEmptyValue("boundary_wall_condition_name").SetString(self.sbm_wall_condition_name)
        settings.AddEmptyValue("extension_operator_type").SetString(self.settings["formulation"]["extension_operator_type"].GetString())
        settings.AddEmptyValue("mls_extension_operator_order").SetInt(self.settings["formulation"]["mls_extension_operator_order"].GetInt())

        # Create a boundary utility for all skin model part names
        time_prev = datetime.datetime.now().replace(microsecond=0)
        self.sbm_utilities = []
        for skin_model_part_name, enclosed_area in zip(self.skin_model_part_names, self.enclosed_areas):
            # Adapt settings
            settings.AddEmptyValue("skin_model_part_name").SetString(skin_model_part_name)
            settings.AddEmptyValue("enclosed_area").SetString(enclosed_area)

            # Create boundary utility
            sbm_utility = KM.ShiftedBoundaryPointBasedUtility(self.model, settings)
            self.sbm_utilities.append(sbm_utility)
            KM.Logger.PrintInfo(self.__class__.__name__, "New shifted-boundary point-based utility created for skin model part '" + skin_model_part_name + "'.")

        if len(self.sbm_utilities) == 0:
            KM.Logger.PrintWarning('No shifted-boundary interface utility was created because no skin model part name was given.')

        self.__PrintAndResetTimer(time_prev, "create the shifted-boundary utilities")

    def __FlagBoundaryElements(self):
        # Flag elements as BOUNDARY if they are located where the skin model parts are
        time_prev = datetime.datetime.now().replace(microsecond=0)
        for sbm_utility in self.sbm_utilities:
            sbm_utility.FindElementsAtTessellatedBoundary()
            sbm_utility.FlagBoundaryElements()
        self.__PrintAndResetTimer(time_prev, "flag boundary elements")

    def __ImmerseShiftedBoundaries(self):
        if len(self.sbm_utilities) > 0:
            # Calculate the required neighbors
            elemental_neighbors_process = KM.GenericFindElementalNeighboursProcess(self.main_model_part)
            elemental_neighbors_process.Execute()

            time_prev = datetime.datetime.now().replace(microsecond=0)

            # Boundary and interface flags must be reset for the volume/ computing model part once before skin model parts are immersed
            self.sbm_utilities[0].ResetFlags()
            time_prev = self.__PrintAndResetTimer(time_prev, "reset flags")

            # Set boundary flags and locate skin model part points in the volume model part elements for all skin model parts.
            # NOTE that boundary elements are flagged after locating the skin points in case skin points are located in elements,
            # which are not touching or intersected by the tessellated boundary
            for sbm_utility in self.sbm_utilities:
                sbm_utility.FindElementsAtTessellatedBoundary()
                time_prev = self.__PrintAndResetTimer(time_prev, "find boundary elements of tessellated skin")
                sbm_utility.MapSkinPointsToElements()
                time_prev = self.__PrintAndResetTimer(time_prev, "find skin points")
                sbm_utility.FlagBoundaryElements()
                time_prev = self.__PrintAndResetTimer(time_prev, "flag boundary elements")

            # Flag interface elements after all boundary elements containing all skin geometries have been found.
            self.sbm_utilities[0].FlagInterfaceElements()
            time_prev = self.__PrintAndResetTimer(time_prev, "flag interface elements")

            # Deactivate BOUNDARY elements and nodes which are surrounded by deactivated elements. Also find and deactivate unstable clusters if requested.
            # An unstable cluster is defined as enclosed fluid volume created by deactivated elements, in which no degree of freedom is fixed.
            #NOTE Right now all cluster except for the biggest one will be deactivated.
            #TODO Separate process for deactivating unstable regions?!
            self.sbm_utilities[0].DeactivateElementsAndNodes(self.deactivate_unstable_clusters)
            time_prev = self.__PrintAndResetTimer(time_prev, "deactivate elements and nodes")

            # Add shifted-boundary conditions for points at the boundary based on extension operators.
            # NOTE that the same boundary sub model part is being used here for all skin model parts and their utilities to add conditions.
            for i_skin, sbm_utility in enumerate(self.sbm_utilities):
                sbm_utility.CalculateAndAddSkinIntegrationPointConditions()
                KM.Logger.PrintInfo(self.__class__.__name__, "Integration point conditions created for skin model part '" + self.skin_model_part_names[i_skin] + "'.")
            time_prev = self.__PrintAndResetTimer(time_prev, "add conditions for all skin model parts")

    def __ReImmerseShiftedBoundaries(self):
        if len(self.sbm_utilities) > 0:
            # Boundary and interface flags must be reset for the volume/ computing model part once before skin model parts are re-immersed
            self.sbm_utilities[0].ResetFlags()

            # Remove the previous conditions from the boundary sub model part
            boundary_sub_model_part = self.GetComputingModelPart().GetSubModelPart(self.boundary_sub_model_part_name)
            for cond in boundary_sub_model_part.Conditions:
                self.main_model_part.GetCondition(cond.Id).Set(KM.TO_ERASE, True)
            self.main_model_part.RemoveConditions(KM.TO_ERASE)

            # Set boundary flags and locate skin model part points in the volume model part elements for all skin model parts.
            for sbm_utility in self.sbm_utilities:
                sbm_utility.UpdateBoundaryElements()
                sbm_utility.MapSkinPointsToElements()
                sbm_utility.FlagBoundaryElements()

            # Flag interface elements after all boundary elements containing all skin geometries have been found.
            self.sbm_utilities[0].FlagInterfaceElements()

            # Deactivate BOUNDARY elements and nodes which are surrounded by deactivated elements. Also find and deactivate unstable clusters if requested.
            # An unstable cluster is defined as enclosed fluid volume created by deactivated elements, in which no degree of freedom is fixed.
            #NOTE Right now all cluster except for the biggest one will be deactivated.
            self.sbm_utilities[0].DeactivateElementsAndNodes(self.deactivate_unstable_clusters)

            # Add shifted-boundary conditions for points at the boundary based on extension operators.
            # NOTE that the same boundary sub model part is being used here for all skin model parts and their utilities to add conditions.
            for sbm_utility in self.sbm_utilities:
                sbm_utility.CalculateAndAddSkinIntegrationPointConditions()

            KM.Logger.PrintInfo(self.__class__.__name__, "Shifted boundaries were re-immersed.")

    def __FillSkinPointsFromDiscModelPart(self, skin_model_part_name, skin_model_part, skin_points_model_part):
        # Create skin points in skin points model part from the discretized skin model part if the skin points model part is empty and the skin model part is not empty.
        n_skin_points = skin_points_model_part.NumberOfNodes()

        # Get integration points of each skin element and add as skin point
        for element in skin_model_part.Elements:
            integration_points = element.GetIntegrationPoints()
            integration_points_weights = element.GetIntegrationPointWeights()
            elem_normal = element.GetGeometry().Normal()
            elem_unit_normal = elem_normal / max(1e-10, numpy.linalg.norm(elem_normal))

            for int_pt, int_pt_w in zip(integration_points, integration_points_weights):
                # Add integration point as node to skin points model part and add its area normal as nodal variable
                new_node = skin_points_model_part.CreateNewNode(n_skin_points, int_pt[0], int_pt[1], int_pt[2])
                n_skin_points += 1

                # Calculate and set the area normal of the integration point
                # NOTE that the length of the normal will be used as integration weight for the boundary condition!
                skin_pt_area_normal = elem_unit_normal * int_pt_w
                new_node.SetValue(KM.NORMAL, skin_pt_area_normal)

        KM.Logger.PrintInfo(self.__class__.__name__, "Skin points created as nodes in skin points model part for skin model part '" + skin_model_part_name + "'.")

    def __UpdateSkinPointsFromDiscModelParts(self):
        # Update skin points in skin points model part from the discretized skin model part if the skin points model part is not empty and the skin model part is not empty.
        # NOTE that this function assumes that the number of skin points in the skin points model part stays equal to the number of integration points in the skin model part.
        for skin_mp_name in self.skin_model_part_names:
            skin_mp = self.model.GetModelPart(skin_mp_name)
            skin_mp_points = self.model.GetModelPart(skin_mp_name+"Points")
            if skin_mp_points.NumberOfNodes() != 0 and skin_mp.NumberOfNodes() != 0:
                # Get integration points of each skin element and update the corresponding node in the skin points model part
                n_skin_points = 0

                for element in skin_mp.Elements:
                    integration_points = element.GetIntegrationPoints()
                    integration_points_weights = element.GetIntegrationPointWeights()
                    elem_normal = element.GetGeometry().Normal()
                    elem_unit_normal = elem_normal / max(1e-10, numpy.linalg.norm(elem_normal))

                    # Get skin velocity at the skin element integration points
                    embedded_velocities = element.CalculateOnIntegrationPoints(KM.VELOCITY, skin_mp.ProcessInfo)

                    for int_pt, int_pt_w, int_pt_u in zip(integration_points, integration_points_weights, embedded_velocities):
                        # Update integration point node coordinates in skin points model part
                        int_pt_node = skin_mp_points.GetNode(n_skin_points)
                        int_pt_node.X = int_pt[0]
                        int_pt_node.Y = int_pt[1]
                        int_pt_node.Z = int_pt[2]

                        # Set the embedded velocity at the integration point
                        int_pt_node.SetValue(KM.EMBEDDED_VELOCITY, int_pt_u)

                        # Calculate and set the area normal of the integration point
                        # NOTE that the length of the normal will be used as integration weight for the boundary condition!
                        skin_pt_area_normal = elem_unit_normal * int_pt_w
                        int_pt_node.SetValue(KM.NORMAL, skin_pt_area_normal)

                        n_skin_points += 1

                KM.Logger.PrintInfo(self.__class__.__name__, "Skin points updated as nodes in skin points model part for skin model part '" + skin_mp_name + "'.")

    def _FmAleIsActive(self):
        return self.settings["fm_ale_settings"]["fm_ale_step_frequency"].GetInt() > 0

    def __GetFmAleVirtualModelPart(self):
        if not hasattr(self, '_virtual_model_part'):
            self._virtual_model_part = self.__CreateFmAleVirtualModelPart()
        return self._virtual_model_part

    def __CreateFmAleVirtualModelPart(self):
        virtual_model_part_name = self.settings["fm_ale_settings"]["fm_ale_solver_settings"]["virtual_model_part_name"].GetString()
        virtual_model_part = self.model.CreateModelPart(virtual_model_part_name)
        return virtual_model_part

    def __GetFmAleUtility(self):
        if not hasattr (self, '_mesh_moving_util'):
            self._mesh_moving_util = self.__CreateFmAleUtility()
        return self._mesh_moving_util

    def __CreateFmAleUtility(self):
        if mesh_moving_available:
            mesh_movement = self.settings["fm_ale_settings"]["mesh_movement"].GetString()
            if (mesh_movement == "implicit"):
                return KM_MeshMoving.FixedMeshALEUtilities(self.model, self.settings["fm_ale_settings"]["fm_ale_solver_settings"])
            else:
                raise Exception("FM-ALE mesh_movement set to \'" + mesh_movement + "\'. Available option is \'implicit\'.")
        else:
            raise Exception("MeshMovingApplication is required to construct the FM-ALE utility (FixedMeshALEUtilities)")

    def __DoFmAleMeshMovementAndProjection(self):
        # Solve the mesh problem
        self.__GetFmAleUtility().ComputeMeshMovement(self.main_model_part.ProcessInfo[KM.DELTA_TIME])

        # Project the obtained MESH_VELOCITY and historical VELOCITY and PRESSURE values to the origin mesh
        buffer_size = self.main_model_part.GetBufferSize()
        domain_size = self.main_model_part.ProcessInfo[KM.DOMAIN_SIZE]
        if (domain_size == 2):
            self.__GetFmAleUtility().ProjectVirtualValues2D(self.main_model_part, buffer_size)
        else:
            self.__GetFmAleUtility().ProjectVirtualValues3D(self.main_model_part, buffer_size)

    def __IsFmAleStep(self):
        if self._FmAleIsActive():
            if (self.fm_ale_step == self.settings["fm_ale_settings"]["fm_ale_step_frequency"].GetInt()):
                return True
            else:
                return False
        else:
            return False

    def __PrintAndResetTimer(self, time_prev, process_description):
        time_curr = datetime.datetime.now().replace(microsecond=0)
        delta = time_curr - time_prev
        hours, remainder = divmod(delta.seconds, 3600)
        minutes, seconds = divmod(remainder, 60)
        KM.Logger.PrintInfo(self.__class__.__name__, "Time required to " + process_description + ": " + str(int(hours)) + " h, " + str(int(minutes)) + " min, " + str(int(seconds)) + " sec")
        return time_curr
