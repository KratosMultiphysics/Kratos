# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.GeoMechanicsApplication as KratosGeo

# Import base class file
from KratosMultiphysics.GeoMechanicsApplication.geomechanics_solver import GeoMechanicalSolver as GeoSolver

def CreateSolver(model, custom_settings):
    return PwSolver(model, custom_settings)

class PwSolver(GeoSolver):
    '''Solver for the solution of pore pressure problems.'''

    def __init__(self, model, custom_settings):
        super().__init__(model, custom_settings)

        KratosMultiphysics.Logger.PrintInfo("GeoMechanics_Pw_Solver", "Construction of Solver finished.")

    @classmethod
    def GetDefaultParameters(cls):
        this_defaults = KratosMultiphysics.Parameters("""{
            "solver_type": "geomechanics_U_Pw_solver",
            "model_part_name": "PorousDomain",
            "domain_size": 2,
            "model_import_settings":{
                "input_type": "mdpa",
                "input_filename": "unknown_name"
            },
            "material_import_settings" :{
                "materials_filename": ""
            },
            "time_stepping": {
                "time_step": 0.1
            },
            "buffer_size": 2,
            "echo_level": 0,
            "reform_dofs_at_each_step": false,
            "clear_storage": false,
            "compute_reactions": false,
            "move_mesh_flag": false,
            "nodal_smoothing": false,
            "solution_type": "quasi_static",
            "scheme_type": "Newmark",
            "newmark_beta": 0.25,
            "newmark_gamma": 0.5,
            "newmark_theta": 0.5,
            "rayleigh_m": 0.0,
            "rayleigh_k": 0.0,
            "strategy_type": "newton_raphson",
            "convergence_criterion": "water_pressure_criterion",
            "water_pressure_relative_tolerance": 1.0e-4,
            "water_pressure_absolute_tolerance": 1.0e-9,
            "residual_relative_tolerance": 1.0e-4,
            "residual_absolute_tolerance": 1.0e-9,
            "desired_iterations"         : 4,
            "max_radius_factor"          : 20.0,
            "min_radius_factor"          : 0.5,
            "max_iterations"             : 15,
            "min_iterations"             : 6,
            "number_cycles"              : 5,
            "increase_factor"            : 2.0,
            "reduction_factor"           : 0.5,
            "calculate_reactions"        : true,
            "max_line_search_iterations" : 5,
            "first_alpha_value"          : 0.5,
            "second_alpha_value"         : 1.0,
            "min_alpha"                  : 0.1,
            "max_alpha"                  : 2.0,
            "line_search_tolerance"      : 0.5,
            "rotation_dofs"              : false,
            "block_builder"              : true,
            "search_neighbours_step"     : false,
            "linear_solver_settings":{
                "solver_type": "amgcl",
                "tolerance": 1.0e-6,
                "max_iteration": 100,
                "scaling": false,
                "verbosity": 0,
                "preconditioner_type": "amg",
                "smoother_type": "ilu0",
                "krylov_type": "gmres",
                "coarsening_type": "aggregation"
            },
            "problem_domain_sub_model_part_list": [""],
            "processes_sub_model_part_list": [""],
            "body_domain_sub_model_part_list": [""],
            "loads_variable_list": []
        }""")

        this_defaults.AddMissingParameters(super().GetDefaultParameters())
        return this_defaults

    def PrepareModelPart(self):
        super().PrepareModelPart()
        KratosMultiphysics.Logger.PrintInfo("GeoMechanics_Pw_Solver", "Model reading finished.")

    def AddDofs(self):
        ## Fluid D.o.F.
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.WATER_PRESSURE, KratosMultiphysics.REACTION_WATER_PRESSURE,self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosGeo.DT_WATER_PRESSURE, self.main_model_part)

        KratosMultiphysics.Logger.PrintInfo("GeoMechanics_Pw_Solver", "DOFs added correctly.")

    def Initialize(self):
        KratosMultiphysics.Logger.PrintInfo("::[GeoMechanics_Pw_Solver]:: ", "Initialisation ...")
        
        super().Initialize()

        KratosMultiphysics.Logger.PrintInfo("GeoMechanics_Pw_Solver", "solver.Initialize is set successfully")

        # Check if everything is assigned correctly
        self.Check()

        KratosMultiphysics.Logger.PrintInfo("GeoMechanics_Pw_Solver", "Solver initialization finished.")

    def _ConstructScheme(self, scheme_type, solution_type):

        self.main_model_part.ProcessInfo.SetValue(KratosGeo.VELOCITY_COEFFICIENT,    1.0)
        self.main_model_part.ProcessInfo.SetValue(KratosGeo.DT_PRESSURE_COEFFICIENT, 1.0)

        if not(solution_type.lower() == "transient-groundwater-flow"    or solution_type.lower() == "transient_groundwater_flow"    or
               solution_type.lower() == "steady-state-groundwater-flow" or solution_type.lower() == "steady_state_groundwater_flow"   ):
            err_msg = "Undefined solution type:" + solution_type + " , only Transient groundwater flow and Steady state groundwater flow are available."
            raise RuntimeError(err_msg)
        KratosMultiphysics.Logger.PrintInfo("GeoMechanics_Pw_Solver, solution_type", solution_type)

        if scheme_type.lower() == "newmark":
            KratosMultiphysics.Logger.PrintInfo("GeoMechanics_Pw_Solver, scheme", "Newmark.")
            return KratosGeo.NewmarkQuasistaticPwScheme(self.settings["newmark_theta"].GetDouble())

        if scheme_type.lower() == "backward_euler":
            KratosMultiphysics.Logger.PrintInfo("GeoMechanics_Pw_Solver, scheme", "Backward Euler.")
            return KratosGeo.BackwardEulerQuasistaticPwScheme()

        raise RuntimeError("Apart from Newmark and Backward Euler, other scheme_type are not available.")

    def _ConstructConvergenceCriterion(self, convergence_criterion):
        if convergence_criterion.lower() == "water_pressure_criterion":
            return self._MakeWaterPressureCriterion()

        if convergence_criterion.lower() == "residual_criterion":
            return self._MakeResidualCriterion()

        if convergence_criterion.lower() == "and_criterion":
            return KratosMultiphysics.AndCriteria(self._MakeResidualCriterion(), self._MakeWaterPressureCriterion())

        if convergence_criterion.lower() == "or_criterion":
            return KratosMultiphysics.OrCriteria(self._MakeResidualCriterion(), self._MakeWaterPressureCriterion())

        err_msg =  "The requested convergence criterion \"" + convergence_criterion + "\" is not available!\n"
        err_msg += "Available options are: \"water_pressure_criterion\", \"residual_criterion\", \"and_criterion\", \"or_criterion\""
        raise RuntimeError(err_msg)
