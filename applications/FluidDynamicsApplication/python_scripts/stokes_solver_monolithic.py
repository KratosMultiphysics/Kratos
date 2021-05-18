# importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.FluidDynamicsApplication

# Import base class file
from KratosMultiphysics.FluidDynamicsApplication.navier_stokes_solver_vmsmonolithic import StabilizedFormulation
from KratosMultiphysics.FluidDynamicsApplication.navier_stokes_solver_vmsmonolithic import NavierStokesSolverMonolithic

class StabilizedStokesFormulation(StabilizedFormulation):
    """
    Helper class to define Stokes formulation parameters.

    Attributes:

    element_name -- The registered element name for the chosen formulation
    condition_neme  -- The registered condition name for the chosen formulation
    element_integrates_in_time -- States if the time integration is wether done by the element or not
    element_has_nodal_properties -- States if the material properties are wether stored in the nodes or not
    process_data -- Auxiliary container to temporary store the elemental ProcessInfo data
    """

    def __init__(self,settings):
        """
        Constructs the Stokes stabilized formulation auxiliary class

        Arguments:

        settings -- Kratos parameter object encapsulating the formulation settings

        """

        self.element_name = None
        self.condition_name = "NavierStokesWallCondition"
        self.element_integrates_in_time = False
        self.element_has_nodal_properties = False
        self.process_data = {}

        if settings.Has("element_type"):
            formulation = settings["element_type"].GetString()
            if formulation == "symbolic_stokes":
                self._SetUpSymbolicStokes(settings)
            else:
                err_msg = "Found \'element_type\' is \'" + formulation + "\'.\n"
                err_msg += "Available options are:\n"
                err_msg += "\t- \'symbolic_stokes\'"
                raise RuntimeError(err_msg)
        else:
            print(settings)
            raise RuntimeError("Argument \'element_type\' not found in stabilization settings.")

    def _SetUpSymbolicStokes(self,settings):
        default_settings = KratosMultiphysics.Parameters(r"""{
            "element_type": "symbolic_stokes",
            "dynamic_tau": 1.0
        }""")
        settings.ValidateAndAssignDefaults(default_settings)

        self.element_name = "SymbolicStokes"
        self.element_integrates_in_time = True

        self.process_data[KratosMultiphysics.DYNAMIC_TAU] = settings["dynamic_tau"].GetDouble()

def CreateSolver(main_model_part, custom_settings):
    return StokesSolverMonolithic(main_model_part, custom_settings)

class StokesSolverMonolithic(NavierStokesSolverMonolithic):
    """
    Monolithic Stokes formulations solver.

    This solver is an specialization of the Navier-Stokes monolithic solver
    designed to be used in combination with elements implementing Stokes formulations.

    The main differences with its base class are the possibility of forcing the steady
    state by deactivating the time integration coefficients and the use of a linear
    strategy by default. The non-linear Newton-Raphson strategy is also supported for
    the non-newtonian material models case.
    """

    @classmethod
    def GetDefaultParameters(cls):
        ##settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
        {
            "solver_type": "stokes_solver_monolithic",
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
                "element_type": "symbolic_stokes"
            },
            "echo_level": 0,
            "compute_reactions": false,
            "analysis_type": "linear",
            "reform_dofs_at_each_step": false,
            "consider_periodic_conditions": false,
            "relative_velocity_tolerance": 1e-3,
            "absolute_velocity_tolerance": 1e-5,
            "relative_pressure_tolerance": 1e-3,
            "absolute_pressure_tolerance": 1e-5,
            "linear_solver_settings"        : {
                "solver_type" : "amgcl"
            },
            "volume_model_part_name" : "volume_model_part",
            "skin_parts": [""],
            "assign_neighbour_elements_to_conditions": true,
            "no_skin_parts":[""],
            "time_stepping"                : {
                "automatic_time_step" : false,
                "CFL_number"          : 1,
                "minimum_delta_time"  : 1e-4,
                "maximum_delta_time"  : 0.01,
                "time_step"           : 0.0
            },
            "time_scheme": "bdf2",
            "force_steady_state": false,
            "move_mesh_flag": false
        }""")

        default_settings.AddMissingParameters(super(StokesSolverMonolithic, cls).GetDefaultParameters())
        return default_settings

    def __init__(self, model, custom_settings):
        """
        Constructs the Stokes monolithic solver

        Arguments:

        model -- The model container to be used
        custom_settings -- Kratos parameter object encapsulating the solver settings

        """

        self._validate_settings_in_baseclass = True # To be removed eventually
        super().__init__(model, custom_settings)
        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Construction of StokesSolverMonolithic finished.")

    def AddVariables(self):
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PRESSURE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION_WATER_PRESSURE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.BODY_FORCE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.EXTERNAL_PRESSURE)

        #TODO: These are required by the NavierStokesWallCondition
        #TODO: ACCELERATION is not required
        #TODO: MESH_VELOCITY is required in case a wall-law is used
        #TODO: Remove as soon as we clean up the conditions
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ACCELERATION)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.MESH_VELOCITY)

        # Adding variables required for the nodal material properties
        if self.element_has_nodal_properties:
            self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DENSITY)
            self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VISCOSITY)

        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Fluid solver variables added correctly.")

    def InitializeSolutionStep(self):
        # Computes BDF coefficients and strategy InitializeSolutionStep
        super().InitializeSolutionStep()

        # Force the steady regime by setting the BDF coefficients to zero
        if(self.settings["force_steady_state"].GetBool()):
            bdf_vec = self.GetComputingModelPart().ProcessInfo[KratosMultiphysics.BDF_COEFFICIENTS]
            for i in range(len(bdf_vec)):
                bdf_vec[i] = 0.0
            self.GetComputingModelPart().ProcessInfo.SetValue(KratosMultiphysics.BDF_COEFFICIENTS, bdf_vec)

    def _SetFormulation(self):
        self.formulation = StabilizedStokesFormulation(self.settings["formulation"])
        self.element_name = self.formulation.element_name
        self.condition_name = self.formulation.condition_name
        self.element_integrates_in_time = self.formulation.element_integrates_in_time
        self.element_has_nodal_properties = self.formulation.element_has_nodal_properties

    def _SetTimeSchemeBufferSize(self):
        scheme_type = self.settings["time_scheme"].GetString()
        if scheme_type == "bdf2":
            self.min_buffer_size = 3
        else:
            msg  = "Unknown time_scheme option found in project parameters:\n"
            msg += "\"" + scheme_type + "\"\n"
            msg += "Accepted value is \"bdf2\".\n"
            raise Exception(msg)
