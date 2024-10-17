# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD

# Import base class file
from KratosMultiphysics.FluidDynamicsApplication.fluid_solver import FluidSolver

class StabilizedFormulation:
    """Helper class to define stabilization-dependent parameters."""
    def __init__(self,settings):
        self.element_name = None
        self.element_integrates_in_time = False
        self.element_has_nodal_properties = False
        self.historical_nodal_properties_variables_list = []
        self.non_historical_nodal_properties_variables_list = []
        self.process_data = {}

        #TODO: Keep this until the MonolithicWallCondition is removed to ensure backwards compatibility in solvers with no defined condition_name
        self.condition_name = "MonolithicWallCondition"

        if settings.Has("element_type"):
            formulation = settings["element_type"].GetString()
            if formulation == "vms":
                self._SetUpClassicVMS(settings)
            elif formulation == "qsvms":
                self._SetUpQSVMS(settings)
            elif formulation == "dvms":
                self._SetUpDVMS(settings)
            elif formulation == "fic":
                self._SetUpFIC(settings)
            elif formulation == "symbolic":
                warn_msg  = 'Provided \'element_name\' is \'symbolic\'. This is has been renamed to \'weakly_compressible\'. Use this instead.'
                KratosMultiphysics.Logger.PrintWarning('\n\x1b[1;31mDEPRECATION-WARNING\x1b[0m', warn_msg)
                self._SetUpWeaklyCompressible(settings)
            elif formulation == "weakly_compressible":
                self._SetUpWeaklyCompressible(settings)
            elif formulation == "axisymmetric_navier_stokes":
                self._SetUpAxisymmetricNavierStokes(settings)
            elif formulation == "p2p1":
                self._SetUpP2P1(settings)
            else:
                formulation_list = ["qsvms", "dvms", "fic", "weakly_compressible", "axisymmetric_navier_stokes", "p2p1"]
                err_msg = f"Wrong \'element_type\' : \'{formulation}\' provided. Available options are:\n"
                for elem in formulation_list:
                    err_msg += f"\t- {elem}\n"
                # raise RuntimeError(err_msg) #TODO: Turn this into an error once the derived solvers handle this properly
                KratosMultiphysics.Logger.PrintWarning("NavierStokesMonolithicSolver", err_msg)
        else:
            print(settings)
            raise RuntimeError("Argument \'element_type\' not found in formulation settings.")

    def SetProcessInfo(self,model_part):
        for variable,value in self.process_data.items():
            model_part.ProcessInfo[variable] = value

    def _SetUpClassicVMS(self,settings):
        default_settings = KratosMultiphysics.Parameters(r"""{
            "element_type": "vms",
            "use_orthogonal_subscales": false,
            "dynamic_tau": 0.01
        }""")

        default_non_newtonian_settings = KratosMultiphysics.Parameters(r"""{
            "power_law_k": 1e-6,
            "power_law_n": 1.0,
            "yield_stress": 0.0,
            "regularization_coefficient" : 100.0
        }""")

        # if non-newtonian, there are some extra options
        if settings.Has("non_newtonian_fluid_parameters"):
            self.non_newtonian_option = True
            default_settings.AddValue("non_newtonian_fluid_parameters", default_non_newtonian_settings)
            self.element_name = 'HerschelBulkleyVMS'
        else:
            self.non_newtonian_option = False
            self.element_name = 'VMS'
        self.condition_name = "MonolithicWallCondition"

        settings.ValidateAndAssignDefaults(default_settings)

        # set the nodal material properties flag
        self.element_has_nodal_properties = True
        self.historical_nodal_properties_variables_list = [KratosMultiphysics.DENSITY, KratosMultiphysics.VISCOSITY]

        # validate the non-newtonian parameters if necessary
        if self.non_newtonian_option:
            settings["non_newtonian_fluid_parameters"].ValidateAndAssignDefaults(default_non_newtonian_settings)
            self.process_data[KratosMultiphysics.POWER_LAW_K] = settings["non_newtonian_fluid_parameters"]["power_law_k"].GetDouble()
            self.process_data[KratosMultiphysics.POWER_LAW_N] = settings["non_newtonian_fluid_parameters"]["power_law_n"].GetDouble()
            self.process_data[KratosMultiphysics.YIELD_STRESS] = settings["non_newtonian_fluid_parameters"]["yield_stress"].GetDouble()
            self.process_data[KratosCFD.REGULARIZATION_COEFFICIENT] = settings["non_newtonian_fluid_parameters"]["regularization_coefficient"].GetDouble()

        self.process_data[KratosMultiphysics.DYNAMIC_TAU] = settings["dynamic_tau"].GetDouble()
        use_oss = settings["use_orthogonal_subscales"].GetBool()
        self.process_data[KratosMultiphysics.OSS_SWITCH] = int(use_oss)

    def _SetUpQSVMS(self,settings):
        default_settings = KratosMultiphysics.Parameters(r"""{
            "element_type": "qsvms",
            "use_orthogonal_subscales": false,
            "dynamic_tau": 0.0,
            "element_manages_time_integration": false
        }""")
        settings.ValidateAndAssignDefaults(default_settings)

        if settings["element_manages_time_integration"].GetBool() == False:
            self.element_name = "QSVMS"
            self.element_integrates_in_time = False
        else:
            self.element_name = "TimeIntegratedQSVMS"
            self.element_integrates_in_time = True
        self.condition_name = "NavierStokesWallCondition"

        self.process_data[KratosMultiphysics.DYNAMIC_TAU] = settings["dynamic_tau"].GetDouble()
        use_oss = settings["use_orthogonal_subscales"].GetBool()
        self.process_data[KratosMultiphysics.OSS_SWITCH] = int(use_oss)

    def _SetUpDVMS(self,settings):
        default_settings = KratosMultiphysics.Parameters(r"""{
            "element_type": "dvms",
            "dynamic_tau": 0.0,
            "use_orthogonal_subscales": false
        }""")
        settings.ValidateAndAssignDefaults(default_settings)

        self.element_name = "DVMS"
        self.condition_name = "NavierStokesWallCondition"

        self.process_data[KratosMultiphysics.DYNAMIC_TAU] = settings["dynamic_tau"].GetDouble()
        use_oss = settings["use_orthogonal_subscales"].GetBool()
        self.process_data[KratosMultiphysics.OSS_SWITCH] = int(use_oss)

    def _SetUpFIC(self,settings):
        default_settings = KratosMultiphysics.Parameters(r"""{
            "element_type": "fic",
            "beta": 0.8,
            "adjust_beta_dynamically": false,
            "dynamic_tau": 0.0
        }""")
        settings.ValidateAndAssignDefaults(default_settings)

        dynamic_beta = settings["adjust_beta_dynamically"].GetBool()
        if dynamic_beta:
            KratosMultiphysics.Logger.PrintWarning("NavierStokesMonolithicSolver","FIC with dynamic beta not yet implemented, using provided beta as a constant value")
        else:
            self.element_name = "FIC"
        self.condition_name = "NavierStokesWallCondition"

        self.process_data[KratosCFD.FIC_BETA] = settings["beta"].GetDouble()
        self.process_data[KratosMultiphysics.OSS_SWITCH] = 0
        self.process_data[KratosMultiphysics.DYNAMIC_TAU] = settings["dynamic_tau"].GetDouble()

    def _SetUpWeaklyCompressible(self,settings):
        #TODO: Remove this after deprecation period is over
        if (settings["element_type"].GetString() == "symbolic"):
            settings["element_type"].SetString("weakly_compressible")

        default_settings = KratosMultiphysics.Parameters(r"""{
            "element_type": "weakly_compressible",
            "dynamic_tau": 1.0,
            "sound_velocity": 1.0e+12
        }""")
        settings.ValidateAndAssignDefaults(default_settings)

        self.element_name = "WeaklyCompressibleNavierStokes"
        self.condition_name = "NavierStokesWallCondition"
        self.element_integrates_in_time = True

        # set the nodal material properties flag
        self.element_has_nodal_properties = True
        self.historical_nodal_properties_variables_list = [KratosMultiphysics.DENSITY]
        self.non_historical_nodal_properties_variables_list = [KratosMultiphysics.SOUND_VELOCITY]

        self.process_data[KratosMultiphysics.DYNAMIC_TAU] = settings["dynamic_tau"].GetDouble()
        #TODO: Remove SOUND_VELOCITY from ProcessInfo. Should be obtained from the properties.
        self.process_data[KratosMultiphysics.SOUND_VELOCITY] = settings["sound_velocity"].GetDouble()

    def _SetUpAxisymmetricNavierStokes(self,settings):
        default_settings = KratosMultiphysics.Parameters(r"""{
            "element_type": "axisymmetric_navier_stokes",
            "dynamic_tau": 1.0
        }""")
        settings.ValidateAndAssignDefaults(default_settings)

        self.element_name = "AxisymmetricNavierStokes"
        self.condition_name = "LineCondition" #TODO: Implement an AxisymmetricNavierStokesCondition
        self.element_integrates_in_time = True

        # set the nodal material properties flag
        self.element_has_nodal_properties = False

        self.process_data[KratosMultiphysics.DYNAMIC_TAU] = settings["dynamic_tau"].GetDouble()

    def _SetUpP2P1(self,settings):
        default_settings = KratosMultiphysics.Parameters(r"""{
            "element_type": "p2p1",
            "dynamic_tau": 1.0
        }""")
        settings.ValidateAndAssignDefaults(default_settings)

        self.element_name = "IncompressibleNavierStokesP2P1Continuous"
        self.condition_name = "NavierStokesP2P1ContinuousWallCondition"
        self.element_integrates_in_time = True
        self.element_has_nodal_properties = False

        self.process_data[KratosMultiphysics.DYNAMIC_TAU] = settings["dynamic_tau"].GetDouble()

def CreateSolver(model, custom_settings):
    return NavierStokesMonolithicSolver(model, custom_settings)

class NavierStokesMonolithicSolver(FluidSolver):

    @classmethod
    def GetDefaultParameters(cls):

        ##settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
        {
            "solver_type": "monolithic",
            "model_part_name": "FluidModelPart",
            "domain_size": -1,
            "model_import_settings": {
                "input_type": "mdpa",
                "input_filename": "unknown_name",
                "reorder": false
            },
            "material_import_settings": {
                "materials_filename": ""
            },
            "formulation": {
                "element_type": "vms"
            },
            "maximum_iterations": 10,
            "echo_level": 0,
            "consider_periodic_conditions": false,
            "builder_and_solver_settings" : {
                "use_block_builder" : true,
                "use_lagrange_BS"   : false,
                "advanced_settings" : { }
            },
            "multi_point_constraints_used": true,
            "compute_reactions": false,
            "analysis_type": "non_linear",
            "reform_dofs_at_each_step": true,
            "assign_neighbour_elements_to_conditions": true,
            "relative_velocity_tolerance": 1e-3,
            "absolute_velocity_tolerance": 1e-5,
            "relative_pressure_tolerance": 1e-3,
            "absolute_pressure_tolerance": 1e-5,
            "linear_solver_settings"        : {
                "solver_type" : "amgcl"
            },
            "volume_model_part_name" : "volume_model_part",
            "skin_parts": [""],
            "no_skin_parts":[""],
            "time_stepping"                : {
                "automatic_time_step" : false,
                "CFL_number"          : 1,
                "minimum_delta_time"  : 1e-4,
                "maximum_delta_time"  : 0.01,
                "time_step"           : 0.0
            },
            "time_scheme":"bossak",
            "alpha":-0.3,
            "velocity_relaxation":0.9,
            "pressure_relaxation":0.9,
            "move_mesh_strategy": 0,
            "periodic": "periodic",
            "move_mesh_flag": false
        }""")

        default_settings.AddMissingParameters(super().GetDefaultParameters())
        return default_settings

    def __init__(self, model, custom_settings):
        super().__init__(model,custom_settings)

        # Set up the auxiliary class with the formulation settings
        self._SetFormulation()

        # Update the default buffer size according to the selected time scheme
        self._SetTimeSchemeBufferSize()

        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Construction of NavierStokesMonolithicSolver finished.")

    def AddVariables(self):
        ## Add base class variables
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ACCELERATION)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.MESH_VELOCITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PRESSURE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.IS_STRUCTURE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.BODY_FORCE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_AREA)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_H)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ADVPROJ)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DIVPROJ)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION_WATER_PRESSURE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.EXTERNAL_PRESSURE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.Y_WALL)
        self.main_model_part.AddNodalSolutionStepVariable(KratosCFD.Q_VALUE)

        # Adding variables required for the nodal material properties
        if self.element_has_nodal_properties:
            for variable in self.historical_nodal_properties_variables_list:
                self.main_model_part.AddNodalSolutionStepVariable(variable)

        # Adding variables required for the periodic conditions
        if self.settings["consider_periodic_conditions"].GetBool() == True:
            self.main_model_part.AddNodalSolutionStepVariable(KratosCFD.PATCH_INDEX)

        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Fluid solver variables added correctly.")

    def Initialize(self):
        # If the solver requires an instance of the stabilized formulation class, set the process info variables
        if hasattr(self, 'formulation'):
            self.formulation.SetProcessInfo(self.GetComputingModelPart())

        # Construct and initialize the solution strategy
        solution_strategy = self._GetSolutionStrategy()
        solution_strategy.SetEchoLevel(self.settings["echo_level"].GetInt())
        solution_strategy.Initialize()

        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Solver initialization finished.")

    def InitializeSolutionStep(self):
        # If required, compute the BDF coefficients
        if hasattr(self, 'time_discretization'):
            (self.time_discretization).ComputeAndSaveBDFCoefficients(self.GetComputingModelPart().ProcessInfo)
        # Perform the solver InitializeSolutionStep
        self._GetSolutionStrategy().InitializeSolutionStep()

    def SolveSolutionStep(self):
        # Call the base fluid solver to solve current time step
        is_converged = super().SolveSolutionStep()

        # If the P2-P1 element is used, postprocess the pressure in the quadratic nodes for the visualization
        # Note that this must be done in here (not in the FinalizeSolutionStep) in case the SolveSolutionStep
        # is called in a non-linear outer loop (e.g. from the FSI or the CHT solvers)
        if self.element_name == "IncompressibleNavierStokesP2P1Continuous":
            KratosCFD.FluidAuxiliaryUtilities.PostprocessP2P1ContinuousPressure(self.GetComputingModelPart())

        return is_converged

    def _SetFormulation(self):
        self.formulation = StabilizedFormulation(self.settings["formulation"])
        self.element_name = self.formulation.element_name
        self.condition_name = self.formulation.condition_name
        self.element_integrates_in_time = self.formulation.element_integrates_in_time
        self.element_has_nodal_properties = self.formulation.element_has_nodal_properties
        self.historical_nodal_properties_variables_list = self.formulation.historical_nodal_properties_variables_list
        self.non_historical_nodal_properties_variables_list = self.formulation.non_historical_nodal_properties_variables_list

    def _SetTimeSchemeBufferSize(self):
        scheme_type = self.settings["time_scheme"].GetString()
        if scheme_type == "bossak":
            self.min_buffer_size = 2
        elif scheme_type == "bdf2":
            self.min_buffer_size = 3
        elif scheme_type == "steady":
            self.min_buffer_size = 1
            self._SetUpSteadySimulation()
        else:
            msg  = "Unknown time_scheme option found in project parameters:\n"
            msg += "\"" + scheme_type + "\"\n"
            msg += "Accepted values are \"bossak\", \"bdf2\" or \"steady\".\n"
            raise Exception(msg)

    def _SetUpSteadySimulation(self):
        '''Overwrite time stepping parameters so that they do not interfere with steady state simulations.'''
        self.settings["time_stepping"]["automatic_time_step"].SetBool(False)
        if self.settings["formulation"].Has("dynamic_tau"):
            self.settings["formulation"]["dynamic_tau"].SetDouble(0.0)

    def _SetNodalProperties(self):
        set_density = KratosMultiphysics.DENSITY in self.historical_nodal_properties_variables_list
        set_viscosity = KratosMultiphysics.VISCOSITY in self.historical_nodal_properties_variables_list
        set_sound_velocity = KratosMultiphysics.SOUND_VELOCITY in self.non_historical_nodal_properties_variables_list

        # Get density and dynamic viscostity from the properties of the first element
        for el in self.main_model_part.Elements:
            # Get DENSITY from properties
            if set_density:
                rho = el.Properties.GetValue(KratosMultiphysics.DENSITY)
                if rho <= 0.0:
                    raise Exception("DENSITY set to {0} in Properties {1}, positive number expected.".format(rho,el.Properties.Id))
            # Get DYNAMIC_VISCOSITY from properties and calculate the kinematic one (VISCOSITY)
            if set_viscosity:
                dyn_viscosity = el.Properties.GetValue(KratosMultiphysics.DYNAMIC_VISCOSITY)
                if dyn_viscosity <= 0.0:
                    raise Exception("DYNAMIC_VISCOSITY set to {0} in Properties {1}, positive number expected.".format(dyn_viscosity,el.Properties.Id))
                kin_viscosity = dyn_viscosity / rho
            # Get SOUND_VELOCITY
            if set_sound_velocity:
                if el.Properties.Has(KratosMultiphysics.SOUND_VELOCITY):
                    sound_velocity = el.Properties.GetValue(KratosMultiphysics.SOUND_VELOCITY)
                else:
                    sound_velocity = 1.0e+12 # Default sound velocity value
                    KratosMultiphysics.Logger.PrintWarning('No \'SOUND_VELOCITY\' value found in Properties {0}. Setting default value {1}'.format(el.Properties.Id, sound_velocity))
                if sound_velocity <= 0.0:
                    raise Exception("SOUND_VELOCITY set to {0} in Properties {1}, positive number expected.".format(sound_velocity, el.Properties.Id))
            break
        else:
            raise Exception("No fluid elements found in the main model part.")

        # Transfer the obtained properties to the nodes
        if set_density:
            KratosMultiphysics.VariableUtils().SetVariable(KratosMultiphysics.DENSITY, rho, self.main_model_part.Nodes)
        if set_viscosity:
            KratosMultiphysics.VariableUtils().SetVariable(KratosMultiphysics.VISCOSITY, kin_viscosity, self.main_model_part.Nodes)
        if set_sound_velocity:
            KratosMultiphysics.VariableUtils().SetNonHistoricalVariable(KratosMultiphysics.SOUND_VELOCITY, sound_velocity, self.main_model_part.Nodes)
