import KratosMultiphysics
import KratosMultiphysics.GeoMechanicsApplication as KratosGeo
from KratosMultiphysics.GeoMechanicsApplication.geomechanics_solver import GeoMechanicalSolver as GeoSolver


def CreateSolver(model, custom_settings):
    return TSolver(model, custom_settings)


class TSolver(GeoSolver):
    """Solver for the solution of thermal problems."""

    def __init__(self, model, custom_settings):
        super().__init__(model, custom_settings)

        KratosMultiphysics.Logger.PrintInfo("GeoMechanics_T_Solver", "Construction of Solver finished.")

    @classmethod
    def GetDefaultParameters(cls):
        this_defaults = KratosMultiphysics.Parameters("""{
            "solver_type": "geomechanics_T_solver",
            "model_part_name": "PorousDomain",
            "domain_size": 2,
            "model_import_settings":{
                "input_type": "mdpa",
                "input_filename": "unknown_name"
            },
            "material_import_settings" :{
                "materials_filename": ""
            },
            "buffer_size": 2,
            "echo_level": 0,
            "rebuild_level": 2,
            "reform_dofs_at_each_step": false,
            "clear_storage": false,
            "compute_reactions": false,
            "move_mesh_flag": false,
            "nodal_smoothing": false,
            "reset_displacements":  false,
            "solution_type": "quasi_static",
            "scheme_type": "Newmark",
            "newmark_beta": 0.25,
            "newmark_gamma": 0.5,
            "newmark_theta": 0.5,
            "strategy_type": "newton_raphson",
            "convergence_criterion": "temperature_criterion",
            "temperature_relative_tolerance": 1.0e-4,
            "temperature_absolute_tolerance": 1.0e-9,
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
                "solver_type": "AMGCL",
                "tolerance": 1.0e-6,
                "max_iteration": 100,
                "scaling": false,
                "verbosity": 0,
                "preconditioner_type": "ILU0Preconditioner",
                "smoother_type": "ilu0",
                "krylov_type": "gmres",
                "coarsening_type": "aggregation"
            },
            "problem_domain_sub_model_part_list": [""],
            "processes_sub_model_part_list": [""],
            "body_domain_sub_model_part_list": [""],
            "loads_sub_model_part_list": [],
            "loads_variable_list": []
        }""")

        this_defaults.AddMissingParameters(super().GetDefaultParameters())
        return this_defaults

    def PrepareModelPart(self):
        super().PrepareModelPart()
        KratosMultiphysics.Logger.PrintInfo("GeoMechanics_T_Solver", "Model reading finished.")

    def AddDofs(self):
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.TEMPERATURE, self.main_model_part)
        KratosMultiphysics.Logger.PrintInfo("GeoMechanics_T_Solver", "DOFs added correctly.")

    def Initialize(self):
        KratosMultiphysics.Logger.PrintInfo("::[GeoMechanics_T_Solver]:: ", "Started Initialization")
        super().Initialize()

        # Check if everything is assigned correctly
        self.Check()

        KratosMultiphysics.Logger.PrintInfo("GeoMechanics_T_Solver", "Solver initialization finished and validated.")

    def _ConstructScheme(self, scheme_type, solution_type):

        self.main_model_part.ProcessInfo.SetValue(KratosGeo.DT_TEMPERATURE_COEFFICIENT, 1.0)

        KratosMultiphysics.Logger.PrintInfo("GeoMechanics_T_Solver, solution_type", solution_type)
        if solution_type.lower() == "transient_heat_transfer":
            if scheme_type.lower() == "newmark" or scheme_type.lower() == "newmark_flow":
                theta = self.settings["newmark_theta"].GetDouble()
                KratosMultiphysics.Logger.PrintInfo("GeoMechanics_T_Solver, scheme", "Newmark Transient heat transfer.")
                scheme = KratosGeo.NewmarkTScheme(theta)
            elif scheme_type.lower() == "backward_euler":
                 KratosMultiphysics.Logger.PrintInfo("GeoMechanics_T_Solver, scheme", "Backward Euler Transient heat transfer.")
                 scheme = KratosGeo.BackwardEulerTScheme()
            else:
                raise RuntimeError("Apart from Newmark and Backward Euler, no other scheme_type is available for thermal calculations.")

        else:
             raise RuntimeError("Undefined solution type", solution_type)

        return scheme



    def _ConstructConvergenceCriterion(self, convergence_criterion):
        if convergence_criterion.lower() == "temperature_criterion":
            return self._MakeTemperatureCriterion()
        elif convergence_criterion.lower() == "residual_criterion":
            return self._MakeResidualCriterion()
        elif convergence_criterion.lower() == "and_criterion":
            temperature = self._MakeTemperatureCriterion()
            residual = self._MakeResidualCriterion()
            return KratosMultiphysics.AndCriteria(residual, temperature)
        elif convergence_criterion.lower() == "or_criterion":
            temperature = self._MakeTemperatureCriterion()
            residual = self._MakeResidualCriterion()
            return KratosMultiphysics.OrCriteria(residual, temperature)
        else:
            err_msg = "The requested convergence criterion \"" + convergence_criterion + "\" is not available!\n"
            err_msg += "Available options are: \"temperature_criterion\", \"residual_criterion\", \"and_criterion\", \"or_criterion\""
            raise RuntimeError(err_msg)

    def _MakeTemperatureCriterion(self):
        relative_tolerance = self.settings["temperature_relative_tolerance"].GetDouble()
        absolute_tolerance = self.settings["temperature_absolute_tolerance"].GetDouble()

        temperature = KratosMultiphysics.MixedGenericCriteria([(KratosMultiphysics.TEMPERATURE, relative_tolerance, absolute_tolerance)])
        temperature.SetEchoLevel(self.settings["echo_level"].GetInt())

        return temperature
