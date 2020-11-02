from __future__ import absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Import applications
from KratosMultiphysics.FluidDynamicsApplication import navier_stokes_solver_vmsmonolithic as NavierMonolithic

# Import base class file
from KratosMultiphysics.SwimmingDEMApplication.fluid_DEM_coupling_solver import FluidDEMSolver

class StabilizedFormulationDEMCoupled(NavierMonolithic.StabilizedFormulation):
    def __init__(self, settings):
        """The default constructor of the class.

        Keyword arguments:
        self -- It signifies an instance of a class.
        settings -- Kratos parameters containing the formulation used (vms, qsvms, etc).
        """
        super(StabilizedFormulationDEMCoupled, self).__init__(settings)
        if settings.Has("element_type"):
            formulation = settings["element_type"].GetString()
            if formulation == "qsvmsDEM":
                self._SetUpQSVMSDEM(settings)
            if formulation == "MonolithicDEM":
                self._SetUpVMSDEM(settings)

    def _SetUpQSVMSDEM(self,settings):
        default_settings = KratosMultiphysics.Parameters(r"""{
            "element_type": "qsvmsDEM",
            "use_orthogonal_subscales": false,
            "dynamic_tau": 0.01,
            "element_name": "QSVMSDEMCoupled"
        }""")
        settings.ValidateAndAssignDefaults(default_settings)

        self.element_name = settings["element_name"].GetString()
        self.element_has_nodal_properties = True

        self.process_data[KratosMultiphysics.DYNAMIC_TAU] = settings["dynamic_tau"].GetDouble()
        use_oss = settings["use_orthogonal_subscales"].GetBool()
        self.process_data[KratosMultiphysics.OSS_SWITCH] = int(use_oss)

    def _SetUpVMSDEM(self,settings):
        default_settings = KratosMultiphysics.Parameters(r"""{
            "element_type": "MonolithicDEM",
            "use_orthogonal_subscales": false,
            "dynamic_tau": 0.01,
            "element_name": "MonolithicDEMCoupled"
        }""")
        settings.ValidateAndAssignDefaults(default_settings)

        self.element_name = settings["element_name"].GetString()
        self.element_has_nodal_properties = True

        self.process_data[KratosMultiphysics.DYNAMIC_TAU] = settings["dynamic_tau"].GetDouble()
        use_oss = settings["use_orthogonal_subscales"].GetBool()
        self.process_data[KratosMultiphysics.OSS_SWITCH] = int(use_oss)

def CreateSolver(model, custom_settings):
    return NavierStokesSolverMonolithicDEM(model, custom_settings)

class NavierStokesSolverMonolithicDEM(FluidDEMSolver, NavierMonolithic.NavierStokesSolverMonolithic):

    @classmethod
    def GetDefaultParameters(cls):

        ##settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
        {
            "solver_type": "MonolithicDEM",
            "model_part_name": "FluidModelPart",
            "domain_size": -1,
            "model_import_settings": {
                "input_type": "mdpa",
                "input_filename": "unknown_name",
                "reorder": false
            },
            "material_import_settings": {
                "materials_filename": "unknown_materials.json"
            },
            "formulation": {
                "element_type": "qsvms"
            },
            "maximum_iterations": 10,
            "echo_level": 0,
            "consider_periodic_conditions": false,
            "compute_reactions": false,
            "analysis_type": "non_linear",
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
            "assign_neighbour_elements_to_conditions": false,
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
            "turbulence_model": "None"
        }""")

        default_settings.AddMissingParameters(super(NavierStokesSolverMonolithicDEM, cls).GetDefaultParameters())
        return default_settings

    def __init__(self, model, custom_settings):
        """The default constructor of the class.

        Keyword arguments:
        self -- It signifies an instance of a class.
        model -- the container of the fluid model part.
        custom_settings -- Kratos parameters containing fluid solver settings.
        """
        self._validate_settings_in_baseclass=True
        custom_settings = self._BackwardsCompatibilityHelper(custom_settings)
        super(NavierStokesSolverMonolithicDEM,self).__init__(model, custom_settings)

        # Set up the auxiliary class with the formulation settings
        self._SetFormulation()

        super(NavierStokesSolverMonolithicDEM,self)._SetTimeSchemeBufferSize()

        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Construction of NavierStokesSolverMonolithicDEM finished.")

    def Initialize(self):
        super(NavierStokesSolverMonolithicDEM, self).Initialize()

        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Solver initialization finished.")

    def InitializeSolutionStep(self):
        super(NavierStokesSolverMonolithicDEM, self).InitializeSolutionStep()

    def _SetFormulation(self):
        self.formulation = StabilizedFormulationDEMCoupled(self.settings["formulation"])
        self.element_name = self.formulation.element_name
        self.condition_name = self.formulation.condition_name
        self.element_integrates_in_time = self.formulation.element_integrates_in_time
        self.element_has_nodal_properties = self.formulation.element_has_nodal_properties