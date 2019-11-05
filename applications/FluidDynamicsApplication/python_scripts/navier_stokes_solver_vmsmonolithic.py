from __future__ import absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD

# Import base class file
from KratosMultiphysics.FluidDynamicsApplication.fluid_solver import FluidSolver

import KratosMultiphysics.python_linear_solver_factory as linear_solver_factory

from KratosMultiphysics.FluidDynamicsApplication.turbulence_model_solver import CreateTurbulenceModel

class StabilizedFormulation(object):
    """Helper class to define stabilization-dependent parameters."""
    def __init__(self,settings):
        self.element_name = None
        self.condition_name = "MonolithicWallCondition"
        self.element_integrates_in_time = False
        self.element_has_nodal_properties = False
        self.process_data = {}

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
                self._SetUpSymbolic(settings)
        else:
            print(settings)
            raise RuntimeError("Argument \'element_type\' not found in stabilization settings.")

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

        settings.ValidateAndAssignDefaults(default_settings)

        # set the nodal material properties flag
        self.element_has_nodal_properties = True

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

        self.process_data[KratosMultiphysics.DYNAMIC_TAU] = settings["dynamic_tau"].GetDouble()
        use_oss = settings["use_orthogonal_subscales"].GetBool()
        self.process_data[KratosMultiphysics.OSS_SWITCH] = int(use_oss)

    def _SetUpDVMS(self,settings):
        default_settings = KratosMultiphysics.Parameters(r"""{
            "element_type": "dvms",
            "use_orthogonal_subscales": false
        }""")
        settings.ValidateAndAssignDefaults(default_settings)

        self.element_name = "DVMS"

        use_oss = settings["use_orthogonal_subscales"].GetBool()
        self.process_data[KratosMultiphysics.OSS_SWITCH] = int(use_oss)


    def _SetUpFIC(self,settings):
        default_settings = KratosMultiphysics.Parameters(r"""{
            "element_type": "fic",
            "beta": 0.8,
            "adjust_beta_dynamically": false
        }""")
        settings.ValidateAndAssignDefaults(default_settings)

        dynamic_beta = settings["adjust_beta_dynamically"].GetBool()
        if dynamic_beta:
            KratosMultiphysics.Logger.PrintWarning("NavierStokesSolverVMSMonolithic","FIC with dynamic beta not yet implemented, using provided beta as a constant value")
        else:
            self.element_name = "FIC"

        self.process_data[KratosCFD.FIC_BETA] = settings["beta"].GetDouble()
        self.process_data[KratosMultiphysics.OSS_SWITCH] = 0

    def _SetUpSymbolic(self,settings):
        default_settings = KratosMultiphysics.Parameters(r"""{
            "element_type": "symbolic",
            "dynamic_tau": 1.0,
            "sound_velocity": 1.0e+12
        }""")
        settings.ValidateAndAssignDefaults(default_settings)

        self.element_name = "SymbolicNavierStokes"
        self.condition_name = "NavierStokesWallCondition"
        self.element_integrates_in_time = True

        self.process_data[KratosMultiphysics.DYNAMIC_TAU] = settings["dynamic_tau"].GetDouble()
        self.process_data[KratosMultiphysics.SOUND_VELOCITY] = settings["sound_velocity"].GetDouble()

def CreateSolver(model, custom_settings):
    return NavierStokesSolverMonolithic(model, custom_settings)

class NavierStokesSolverMonolithic(FluidSolver):

    @classmethod
    def GetDefaultSettings(cls):

        ##settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
        {
            "solver_type": "navier_stokes_solver_vmsmonolithic",
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
            "compute_reactions": false,
            "reform_dofs_at_each_step": true,
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
            "move_mesh_flag": false,
            "turbulence_model_solver_settings": {}
        }""")

        default_settings.AddMissingParameters(super(NavierStokesSolverMonolithic, cls).GetDefaultSettings())
        return default_settings

    def _BackwardsCompatibilityHelper(self,settings):
        ## Backwards compatibility -- deprecation warnings
        if settings.Has("stabilization"):
            msg  = "Input JSON data contains deprecated setting \'stabilization\'.\n"
            msg += "Please rename it to \'formulation\' (and rename \'stabilization/formulation\' to \'formulation/element_type\' if it exists).\n"
            KratosMultiphysics.Logger.PrintWarning("NavierStokesVMSMonolithicSolver",msg)
            settings.AddValue("formulation", settings["stabilization"])
            settings.RemoveValue("stabilization")
            settings["formulation"].AddValue("element_type", settings["formulation"]["formulation"])
            settings["formulation"].RemoveValue("formulation")

        if settings.Has("oss_switch"):
            msg  = "Input JSON data contains deprecated setting \'oss_switch\' (int).\n"
            msg += "Please define \'formulation/element_type\' (set it to \'vms\')\n"
            msg += "and set \'formulation/use_orthogonal_subscales\' (bool) instead."
            KratosMultiphysics.Logger.PrintWarning("NavierStokesVMSMonolithicSolver",msg)
            if not settings.Has("formulation"):
                settings.AddValue("formulation",KratosMultiphysics.Parameters(r'{"element_type":"vms"}'))
            settings["formulation"].AddEmptyValue("use_orthogonal_subscales")
            settings["formulation"]["use_orthogonal_subscales"].SetBool(bool(settings["oss_switch"].GetInt()))
            settings.RemoveValue("oss_switch")

        if settings.Has("dynamic_tau"):
            msg  = "Input JSON data contains deprecated setting \'dynamic_tau\' (float).\n"
            msg += "Please define \'formulation/element_type\' (set it to \'vms\') and \n"
            msg += "set \'formulation/dynamic_tau\' (float) instead."
            KratosMultiphysics.Logger.PrintWarning("NavierStokesVMSMonolithicSolver",msg)
            if not settings.Has("formulation"):
                settings.AddValue("formulation",KratosMultiphysics.Parameters(r'{"element_type":"vms"}'))
            settings["formulation"].AddEmptyValue("dynamic_tau")
            settings["formulation"]["dynamic_tau"].SetDouble(settings["dynamic_tau"].GetDouble())
            settings.RemoveValue("dynamic_tau")

        if settings.Has("turbulence_model") and settings["turbulence_model"].IsString():
            if settings["turbulence_model"].GetString().lower()!="none":
                msg = "Ignoring deprecated \"turbulence_model\" (string) setting."
                KratosMultiphysics.Logger.PrintWarning("NavierStokesVMSMonolithicSolver",msg)
            settings.RemoveValue("turbulence_model")

        return settings


    def __init__(self, model, custom_settings):
        self._validate_settings_in_baseclass=True # To be removed eventually
        custom_settings = self._BackwardsCompatibilityHelper(custom_settings)
        super(NavierStokesSolverMonolithic,self).__init__(model,custom_settings)

        self.formulation = StabilizedFormulation(self.settings["formulation"])
        self.element_name = self.formulation.element_name
        self.condition_name = self.formulation.condition_name
        self.element_integrates_in_time = self.formulation.element_integrates_in_time
        self.element_has_nodal_properties = self.formulation.element_has_nodal_properties

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

        ## Construct the linear solver
        self.linear_solver = linear_solver_factory.ConstructSolver(self.settings["linear_solver_settings"])

        ## Construct the turbulence model solver
        if not self.settings["turbulence_model_solver_settings"].IsEquivalentTo(KratosMultiphysics.Parameters("{}")):
            self._turbulence_model_solver = CreateTurbulenceModel(model, self.settings["turbulence_model_solver_settings"])
            self.condition_name = self._turbulence_model_solver.GetFluidVelocityPressureConditionName()
            KratosMultiphysics.Logger.PrintInfo("NavierStokesSolverMonolithic", "Using " + self.condition_name + " as wall condition")

        KratosMultiphysics.Logger.PrintInfo("NavierStokesSolverMonolithic", "Construction of NavierStokesSolverMonolithic finished.")


    def AddVariables(self):
        ## Add base class variables
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ACCELERATION)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.MESH_VELOCITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PRESSURE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.IS_STRUCTURE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VISCOSITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DENSITY)
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

        # Adding variables required for the turbulence modelling
        if hasattr(self, "_turbulence_model_solver"):
            self._turbulence_model_solver.fluid_model_part = self.main_model_part
            self._turbulence_model_solver.AddVariables()

        if self.settings["consider_periodic_conditions"].GetBool() == True:
            self.main_model_part.AddNodalSolutionStepVariable(KratosCFD.PATCH_INDEX)

        KratosMultiphysics.Logger.PrintInfo("NavierStokesSolverMonolithic", "Fluid solver variables added correctly.")

    def AddDofs(self):
        super(NavierStokesSolverMonolithic, self).AddDofs()

        # Adding DOFs required for the turbulence modelling
        if hasattr(self, "_turbulence_model_solver"):
            self._turbulence_model_solver.AddDofs()

    def PrepareModelPart(self):
        super(NavierStokesSolverMonolithic, self).PrepareModelPart()

        # Missing prepare model part operations required for the turbulence modelling
        if hasattr(self, "_turbulence_model_solver"):
            self._turbulence_model_solver.PrepareModelPart()

    def Initialize(self):

        self.computing_model_part = self.GetComputingModelPart()

        # If needed, create the estimate time step utility
        if (self.settings["time_stepping"]["automatic_time_step"].GetBool()):
            self.EstimateDeltaTimeUtility = self._GetAutomaticTimeSteppingUtility()

        # Creating the solution strategy

        if self.settings["time_scheme"].GetString() == "steady":
            self.conv_criteria = KratosMultiphysics.ResidualCriteria(self.settings["relative_velocity_tolerance"].GetDouble(),
                                                                     self.settings["absolute_velocity_tolerance"].GetDouble())
        else:
            self.conv_criteria = KratosCFD.VelPrCriteria(self.settings["relative_velocity_tolerance"].GetDouble(),
                                                         self.settings["absolute_velocity_tolerance"].GetDouble(),
                                                         self.settings["relative_pressure_tolerance"].GetDouble(),
                                                         self.settings["absolute_pressure_tolerance"].GetDouble())

        (self.conv_criteria).SetEchoLevel(self.settings["echo_level"].GetInt())

        # Creating the time integration scheme
        if (self.element_integrates_in_time):
            # "Fake" scheme for those cases in where the element manages the time integration
            # It is required to perform the nodal update once the current time step is solved
            self.time_scheme = KratosMultiphysics.ResidualBasedIncrementalUpdateStaticSchemeSlip(
                self.computing_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE],
                self.computing_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]+1)
            # In case the BDF2 scheme is used inside the element, the BDF time discretization utility is required to update the BDF coefficients
            if (self.settings["time_scheme"].GetString() == "bdf2"):
                time_order = 2
                self.time_discretization = KratosMultiphysics.TimeDiscretization.BDF(time_order)
            else:
                err_msg = "Requested elemental time scheme \"" + self.settings["time_scheme"].GetString()+ "\" is not available.\n"
                err_msg += "Available options are: \"bdf2\""
                raise Exception(err_msg)
        else:
            if not hasattr(self, "_turbulence_model_solver"):
                # Bossak time integration scheme
                if self.settings["time_scheme"].GetString() == "bossak":
                    if self.settings["consider_periodic_conditions"].GetBool() == True:
                        self.time_scheme = KratosCFD.ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulent(
                            self.settings["alpha"].GetDouble(),
                            self.computing_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE],
                            KratosCFD.PATCH_INDEX)
                    else:
                        self.time_scheme = KratosCFD.ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulent(
                            self.settings["alpha"].GetDouble(),
                            self.settings["move_mesh_strategy"].GetInt(),
                            self.computing_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE])
                # BDF2 time integration scheme
                elif self.settings["time_scheme"].GetString() == "bdf2":
                    self.time_scheme = KratosCFD.GearScheme()
                # Time scheme for steady state fluid solver
                elif self.settings["time_scheme"].GetString() == "steady":
                    self.time_scheme = KratosCFD.ResidualBasedSimpleSteadyScheme(
                            self.settings["velocity_relaxation"].GetDouble(),
                            self.settings["pressure_relaxation"].GetDouble(),
                            self.computing_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE])
                else:
                    err_msg = "Requested time scheme " + self.settings["time_scheme"].GetString() + " is not available.\n"
                    err_msg += "Available options are: \"bossak\", \"bdf2\" and \"steady\""
                    raise Exception(err_msg)
            else:
                self._turbulence_model_solver.Initialize()
                if self.settings["time_scheme"].GetString() == "bossak":
                    self.time_scheme = KratosCFD.ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulent(
                                self.settings["alpha"].GetDouble(),
                                self.settings["move_mesh_strategy"].GetInt(),
                                self.computing_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE],
                                self.settings["turbulence_model_solver_settings"]["velocity_pressure_relaxation_factor"].GetDouble(),
                                self._turbulence_model_solver.GetTurbulenceSolvingProcess())
                # Time scheme for steady state fluid solver
                elif self.settings["time_scheme"].GetString() == "steady":
                    self.time_scheme = KratosCFD.ResidualBasedSimpleSteadyScheme(
                            self.settings["velocity_relaxation"].GetDouble(),
                            self.settings["pressure_relaxation"].GetDouble(),
                            self.computing_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE],
                            self._turbulence_model_solver.GetTurbulenceSolvingProcess())
        if self.settings["consider_periodic_conditions"].GetBool():
            builder_and_solver = KratosCFD.ResidualBasedBlockBuilderAndSolverPeriodic(self.linear_solver,
                                                                                KratosCFD.PATCH_INDEX)
        else:
            builder_and_solver = KratosMultiphysics.ResidualBasedBlockBuilderAndSolver(self.linear_solver)


        self.solver = KratosMultiphysics.ResidualBasedNewtonRaphsonStrategy(self.computing_model_part,
                                                                            self.time_scheme,
                                                                            self.linear_solver,
                                                                            self.conv_criteria,
                                                                            builder_and_solver,
                                                                            self.settings["maximum_iterations"].GetInt(),
                                                                            self.settings["compute_reactions"].GetBool(),
                                                                            self.settings["reform_dofs_at_each_step"].GetBool(),
                                                                            self.settings["move_mesh_flag"].GetBool())

        (self.solver).SetEchoLevel(self.settings["echo_level"].GetInt())

        self.formulation.SetProcessInfo(self.computing_model_part)

        (self.solver).Initialize()

        KratosMultiphysics.Logger.PrintInfo("NavierStokesSolverMonolithic", "Solver initialization finished.")

    def InitializeSolutionStep(self):
        if self._TimeBufferIsInitialized():
            # If required, compute the BDF coefficients
            if hasattr(self, 'time_discretization'):
                (self.time_discretization).ComputeAndSaveBDFCoefficients(self.GetComputingModelPart().ProcessInfo)
            # Perform the solver InitializeSolutionStep
            (self.solver).InitializeSolutionStep()
            # Perform the turbulence modelling InitializeSolutionStep
            if hasattr(self, "_turbulence_model_solver"):
                self._turbulence_model_solver.InitializeSolutionStep()

    def FinalizeSolutionStep(self):
        super(NavierStokesSolverMonolithic, self).FinalizeSolutionStep()
        # Perform the turbulence modelling FinalizeSolutionStep
        if hasattr(self, "_turbulence_model_solver"):
            self._turbulence_model_solver.FinalizeSolutionStep()

    def Check(self):
        super(NavierStokesSolverMonolithic, self).Check()
        # Turbulence modelling check operations
        if hasattr(self, "_turbulence_model_solver"):
            self._turbulence_model_solver.Check()

    def _SetUpSteadySimulation(self):
        '''Overwrite time stepping parameters so that they do not interfere with steady state simulations.'''
        self.settings["time_stepping"]["automatic_time_step"].SetBool(False)
        if self.settings["formulation"].Has("dynamic_tau"):
            self.settings["formulation"]["dynamic_tau"].SetDouble(0.0)