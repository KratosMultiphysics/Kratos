# Importing the Kratos Library
import KratosMultiphysics
from KratosMultiphysics.kratos_utilities import IssueDeprecationWarning

# Import applications
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD

# Import base class file
from KratosMultiphysics.FluidDynamicsApplication.adjoint_fluid_solver import AdjointFluidSolver

class StabilizedAdjointFormulation(object):
    """Helper class to define stabilization-dependent parameters."""
    def __init__(self,settings):
        self.element_name = None
        self.element_integrates_in_time = False
        self.element_has_nodal_properties = False
        self.process_data = {}

        if settings.Has("element_type"):
            formulation = settings["element_type"].GetString()
            if formulation == "vms":
                self._SetUpClassicAdjointVMS(settings)
            elif formulation == "qsvms":
                self._SetUpAdjointQSVMS(settings)
        else:
            print(settings)
            raise RuntimeError("Argument \'element_type\' not found in stabilization settings.")

    def SetProcessInfo(self,model_part):
        for variable,value in self.process_data.items():
            model_part.ProcessInfo[variable] = value

    def _SetUpClassicAdjointVMS(self,settings):
        IssueDeprecationWarning('StabilizedAdjointFormulation', 'Please use the "qsvms" formulation instead of "vms" formulation')
        default_settings = KratosMultiphysics.Parameters(r"""{
            "element_type": "vms",
            "use_orthogonal_subscales": false,
            "dynamic_tau" : 0.0
        }""")

        # if non-newtonian, there are some extra options
        if settings.Has("non_newtonian_fluid_parameters"):
            raise Exception("Adjoints are not yet supported for HerschelBulkleyVMS")

        self.element_name = 'VMSAdjointElement'

        settings.ValidateAndAssignDefaults(default_settings)

        # set the nodal material properties flag
        self.element_has_nodal_properties = True

        self.process_data[KratosMultiphysics.DYNAMIC_TAU] = settings["dynamic_tau"].GetDouble()
        use_oss = settings["use_orthogonal_subscales"].GetBool()
        self.process_data[KratosMultiphysics.OSS_SWITCH] = int(use_oss)

    def _SetUpAdjointQSVMS(self,settings):
        default_settings = KratosMultiphysics.Parameters(r"""{
            "element_type": "qsvms",
            "use_orthogonal_subscales": false,
            "dynamic_tau": 0.0,
            "element_manages_time_integration": false
        }""")
        settings.ValidateAndAssignDefaults(default_settings)

        if settings["element_manages_time_integration"].GetBool() == False:
            self.element_name = "QSVMSAdjoint"
            self.element_integrates_in_time = False
        else:
            raise Exception("Only bossak time integrations is supported with adjoints")

        self.process_data[KratosMultiphysics.DYNAMIC_TAU] = settings["dynamic_tau"].GetDouble()
        use_oss = settings["use_orthogonal_subscales"].GetBool()
        self.process_data[KratosMultiphysics.OSS_SWITCH] = int(use_oss)

def CreateSolver(main_model_part, custom_settings):
    return AdjointMonolithicSolver(main_model_part, custom_settings)

class AdjointMonolithicSolver(AdjointFluidSolver):
    @classmethod
    def GetDefaultParameters(cls):
        # default settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
        {
            "solver_type" : "adjoint_monolithic_solver",
            "model_part_name": "",
            "domain_size": -1,
            "scheme_settings" : {
                "scheme_type" : "bossak"
            },
            "response_function_settings" : {
                "response_type" : "drag"
            },
            "sensitivity_settings" : {},
            "model_import_settings" : {
                "input_type"     : "mdpa",
                "input_filename" : "unknown_name"
            },
            "material_import_settings": {
                "materials_filename": ""
            },
            "linear_solver_settings" : {
                "solver_type" : "amgcl"
            },
            "volume_model_part_name" : "volume_model_part",
            "skin_parts"  : [""],
            "no_skin_parts"  : [""],
            "echo_level"  : 0,
            "time_stepping"               : {
                "automatic_time_step" : false,
                "time_step"           : -0.1
            },
            "consider_periodic_conditions": false,
            "assign_neighbour_elements_to_conditions": true,
            "formulation": {
                "element_type": "qsvms"
            }
        }""")

        default_settings.AddMissingParameters(super().GetDefaultParameters())
        return default_settings

    def __init__(self, model, custom_settings):
        super().__init__(model,custom_settings)
        self.element_has_nodal_properties = True

        self._SetFormulation()

        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Construction of AdjointMonolithicSolver finished.")

    def AddVariables(self):
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.MESH_VELOCITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ACCELERATION)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PRESSURE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL)
        self.main_model_part.AddNodalSolutionStepVariable(KratosCFD.ADJOINT_FLUID_VECTOR_1)
        self.main_model_part.AddNodalSolutionStepVariable(KratosCFD.ADJOINT_FLUID_VECTOR_2)
        self.main_model_part.AddNodalSolutionStepVariable(KratosCFD.ADJOINT_FLUID_VECTOR_3)
        self.main_model_part.AddNodalSolutionStepVariable(KratosCFD.AUX_ADJOINT_FLUID_VECTOR_1)
        self.main_model_part.AddNodalSolutionStepVariable(KratosCFD.ADJOINT_FLUID_SCALAR_1)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VISCOSITY)
        # TODO: Remove when old VMS elements are removed.
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ADVPROJ)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DIVPROJ)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DENSITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.BODY_FORCE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.EXTERNAL_PRESSURE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.SHAPE_SENSITIVITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL_SENSITIVITY)

        if self.settings["consider_periodic_conditions"].GetBool() == True:
            self.main_model_part.AddNodalSolutionStepVariable(KratosCFD.PATCH_INDEX)

        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Adjoint fluid solver variables added correctly.")

    def Initialize(self):
        # If the solver requires an instance of the stabilized formulation class, set the process info variables
        if hasattr(self, 'formulation'):
            self.formulation.SetProcessInfo(self.GetComputingModelPart())

        # Construct and set the solution strategy
        solution_strategy = self._GetSolutionStrategy()
        solution_strategy.SetEchoLevel(self.settings["echo_level"].GetInt())

        # clear RELAXED_ACCELERATION if it is a steady adjoint problem
        scheme_type = self.settings["scheme_settings"]["scheme_type"].GetString()
        if (scheme_type == "steady"):
            KratosMultiphysics.VariableUtils().SetNonHistoricalVariableToZero(KratosMultiphysics.RELAXED_ACCELERATION, self.main_model_part.Nodes)

        # Initialize the strategy and adjoint utilities
        solution_strategy.Initialize()
        self.GetResponseFunction().Initialize()
        self.GetSensitivityBuilder().Initialize()

        domain_size = self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
        KratosMultiphysics.NormalCalculationUtils().CalculateOnSimplex(self.main_model_part.Conditions, domain_size)
        KratosMultiphysics.NormalCalculationUtils().CalculateNormalShapeDerivativesOnSimplex(self.main_model_part.Conditions, domain_size)

        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Solver initialization finished.")

    def _CreateScheme(self):
        response_function = self.GetResponseFunction()
        scheme_type = self.settings["scheme_settings"]["scheme_type"].GetString()
        domain_size = self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]

        # the schemes used in fluid supports SLIP conditions which rotates element/condition
        # matrices based on nodal NORMAL. Hence, the consistent adjoints also need to
        # rotate adjoint element/condition matrices accordingly. Therefore, following
        # schemes are used.
        if scheme_type == "bossak":
            scheme = KratosCFD.VelocityBossakAdjointScheme(self.settings["scheme_settings"], response_function, domain_size, domain_size + 1)
        elif scheme_type == "steady":
            scheme = KratosCFD.SimpleSteadyAdjointScheme(response_function, domain_size, domain_size + 1)
        else:
            raise Exception("Invalid scheme_type: " + scheme_type)
        return scheme

    def _CreateSolutionStrategy(self):
        computing_model_part = self.GetComputingModelPart()
        time_scheme = self._GetScheme()
        builder_and_solver = self._GetBuilderAndSolver()
        calculate_reaction_flag = False
        reform_dof_set_at_each_step = False
        calculate_norm_dx_flag = False
        move_mesh_flag = False
        return KratosMultiphysics.ResidualBasedLinearStrategy(
            computing_model_part,
            time_scheme,
            builder_and_solver,
            calculate_reaction_flag,
            reform_dof_set_at_each_step,
            calculate_norm_dx_flag,
            move_mesh_flag)

    def _SetFormulation(self):
        self.formulation = StabilizedAdjointFormulation(self.settings["formulation"])
        self.element_name = self.formulation.element_name
        self.element_has_nodal_properties = self.formulation.element_has_nodal_properties

        # TODO: at the moment, qs_vms with SLIP adjoints are not supported, only no slip is supported
        #       because AdjointMonolithicWallCondition is based on nodal nu, where as QSVMS is based on
        #       CLs. There need to update navier_stokes_wall_condition.h with adjoints
        self.condition_name = "AdjointMonolithicWallCondition"
