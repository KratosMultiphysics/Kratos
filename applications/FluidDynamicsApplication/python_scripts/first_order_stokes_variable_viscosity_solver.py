# importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.FluidDynamicsApplication as FDApp

# Import base class file
from KratosMultiphysics.FluidDynamicsApplication.navier_stokes_monolithic_solver import StabilizedFormulation
from KratosMultiphysics.FluidDynamicsApplication.navier_stokes_monolithic_solver import NavierStokesMonolithicSolver

class FirstOrderStokesVariableViscosityFormulation(StabilizedFormulation):
    """
    Helper class to define custom Stokes formulation to use in FirstOrderStokesVariableViscositySolver.

    Attributes:

    element_name -- The registered element name for the chosen formulation
    condition_name  -- The registered condition name for the chosen formulation
    element_integrates_in_time -- States if the time integration is whether done by the element or not
    element_has_nodal_properties -- States if the material properties are whether stored in the nodes or not
    process_data -- Auxiliary container to temporary store the elemental ProcessInfo data
    get_conition_neighbours -- Flag to indicate if knowing the conditions neighbouring each node is necessary 
    """

    BVS_GL_ELEMENT = "first_order_stokes_variable_viscosity_bvs_gl"
    PSPG_SD_ELEMENT = "first_order_stokes_variable_viscosity_pspg_sd"

    ELEMENTS_DICT = {
        BVS_GL_ELEMENT: {
            "element_name": "FirstOrderStokesVariableViscosityBvsGl",
            "get_condition_neighbours": True
        },
        PSPG_SD_ELEMENT: {
            "element_name": "FirstOrderStokesVariableViscosityPspgSd",
            "get_condition_neighbours": False
        }
    }
    

    def __init__(self,settings):
        """
        Constructs the Stokes stabilized formulation auxiliary class

        Arguments:

        settings -- Kratos parameter object encapsulating the formulation settings

        """

        self.element_name = None
        self.condition_name = "FirstOrderStokesVariableViscosityCondition"
        self.delta_parameter = None
        self.sigma_parameter = None
        self.element_integrates_in_time = False
        self.element_has_nodal_properties = True
        self.process_data = {}
        self.get_neighbour_conditions = False

        allowed_element_types = self.ELEMENTS_DICT.keys()

        if settings.Has("element_type"):
            element_type = settings["element_type"].GetString()
            if element_type in allowed_element_types:
                self._SetUpCustomStokesVariableViscosityElement(settings)
            else:
                err_msg = "Found \'element_type\' is \'" + element_type + "\'.\n"
                err_msg += "Available options are:\n\t - \'"
                err_msg += "\'\n\t- \'".join(allowed_element_types)
                err_msg += "\'"
                raise RuntimeError(err_msg)
        else:
            print(settings)
            raise RuntimeError("Argument \'element_type\' not found in stabilization settings.")

    def _SetUpCustomStokesVariableViscosityElement(self,settings):
        default_settings = KratosMultiphysics.Parameters(r"""{
            "element_type": "first_order_stokes_variable_viscosity_bvs_gl",
            "delta_parameter": 0,                      
            "sigma_parameter": 0
        }""")
        settings.ValidateAndAssignDefaults(default_settings)

        element_type = settings["element_type"].GetString()
        self.delta_parameter = settings["delta_parameter"].GetDouble()
        self.sigma_parameter = settings["sigma_parameter"].GetDouble()

        self.element_name = self.ELEMENTS_DICT[element_type]["element_name"]
        self.get_neighbour_conditions = self.ELEMENTS_DICT[element_type].get("get_condition_neighbours",False)

def CreateSolver(main_model_part, custom_settings):
    return FirstOrderStokesVariableViscositySolver(main_model_part, custom_settings)

class FirstOrderStokesVariableViscositySolver(NavierStokesMonolithicSolver):
    """
    Stokes Solver to be used with the elements:
        - first_order_stokes_variable_viscosity_bvs_gl

    For more details, check the corresponding directory in the automatic_differentiation folder.
    """

    @classmethod
    def GetDefaultParameters(cls):
        ##settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
        {
            "solver_type": "first_order_stokes_variable_viscosity",
            "model_part_name": "FluidModelPart",
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
            "formulation": {
                "element_type": "first_order_stokes_variable_viscosity_bvs_gl",
                "delta_parameter": 0,
                "sigma_parameter": 0
            },
            "echo_level": 0,
            "compute_reactions": true,
            "analysis_type": "non_linear",
            "reform_dofs_at_each_step": false,
            "consider_periodic_conditions": false,
            "relative_velocity_tolerance": 1e-3,
            "absolute_velocity_tolerance": 1e-5,
            "relative_pressure_tolerance": 1e-3,
            "absolute_pressure_tolerance": 1e-5,
            "linear_solver_settings"        : {
                "solver_type": "LinearSolversApplication.sparse_lu"
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
            "time_scheme": "steady",
            "force_steady_state": true,
            "move_mesh_flag": false
        }""")

        default_settings.AddMissingParameters(super(FirstOrderStokesVariableViscositySolver, cls).GetDefaultParameters())
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
        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Construction of FirstOrderStokesVariableViscositySolver finished.")

    def AddVariables(self):
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PRESSURE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.BODY_FORCE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION_WATER_PRESSURE)

        # Adding variables required for the nodal material properties
        if self.element_has_nodal_properties:
            self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DYNAMIC_VISCOSITY)

        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Fluid solver variables added correctly.")

    def Initialize(self):

        if self.get_neighbour_conditions:
            # Call the FindConditionsNeighours so that, within C++, we might get the conditions associated to each node.
            domain_size = self.settings["domain_size"].GetInt()
            proc = KratosMultiphysics.FindConditionsNeighboursProcess(self.main_model_part, domain_size, 10)
            proc.Execute()

        super().Initialize()

        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "FirstOrderStokesVariableViscositySolver initialization finished.")

    def InitializeSolutionStep(self):
        # Computes BDF coefficients and set delta parameter
        super().InitializeSolutionStep()

        self.GetComputingModelPart().ProcessInfo.SetValue(FDApp.TAUONE, self.delta_parameter)
        self.GetComputingModelPart().ProcessInfo.SetValue(FDApp.TAUTWO, self.sigma_parameter)

        # Force the steady regime by setting the BDF coefficients to zero
        if(self.settings["force_steady_state"].GetBool()):
            bdf_vec = self.GetComputingModelPart().ProcessInfo[KratosMultiphysics.BDF_COEFFICIENTS]
            for i in range(len(bdf_vec)):
                bdf_vec[i] = 0.0
            self.GetComputingModelPart().ProcessInfo.SetValue(KratosMultiphysics.BDF_COEFFICIENTS, bdf_vec)

    def _SetFormulation(self):
        self.formulation = FirstOrderStokesVariableViscosityFormulation(self.settings["formulation"])
        self.element_name = self.formulation.element_name
        self.condition_name = self.formulation.condition_name
        self.delta_parameter = self.formulation.delta_parameter
        self.sigma_parameter = self.formulation.sigma_parameter
        self.element_integrates_in_time = self.formulation.element_integrates_in_time
        self.element_has_nodal_properties = self.formulation.element_has_nodal_properties
        self.get_neighbour_conditions = self.formulation.get_neighbour_conditions

    def _SetNodalProperties(self):
        # Get density and dynamic viscosity from the properties of the first element
        for el in self.main_model_part.Elements:
            # Get DYNAMIC_VISCOSITY from properties and calculate the kinematic one (VISCOSITY)
            dyn_viscosity = el.Properties.GetValue(KratosMultiphysics.DYNAMIC_VISCOSITY)
            if dyn_viscosity <= 0.0:
                raise Exception("DYNAMIC_VISCOSITY set to {0} in Properties {1}, positive number expected.".format(dyn_viscosity,el.Properties.Id))
            break
        else:
            raise Exception("No fluid elements found in the main model part.")

        KratosMultiphysics.VariableUtils().SetVariable(KratosMultiphysics.DYNAMIC_VISCOSITY, dyn_viscosity, self.main_model_part.Nodes)
