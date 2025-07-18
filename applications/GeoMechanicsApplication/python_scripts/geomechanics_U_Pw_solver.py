# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.GeoMechanicsApplication as KratosGeo
import KratosMultiphysics.StructuralMechanicsApplication as KratosStructure

# Import base class file
from KratosMultiphysics.GeoMechanicsApplication.geomechanics_solver import GeoMechanicalSolver as GeoSolver

def CreateSolver(model, custom_settings):
    return UPwSolver(model, custom_settings)

class UPwSolver(GeoSolver):
    '''Solver for the solution of displacement-pore pressure coupled problems.'''

    def __init__(self, model, custom_settings):
        super().__init__(model, custom_settings)

        KratosMultiphysics.Logger.PrintInfo("GeoMechanics_U_Pw_Solver", "Construction of Solver finished.")

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
            "rebuild_level": 2,
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
            "max_piping_iterations": 50,
            "convergence_criterion": "Displacement_criterion",
            "water_pressure_relative_tolerance": 1.0e-4,
            "water_pressure_absolute_tolerance": 1.0e-9,
            "displacement_relative_tolerance": 1.0e-4,
            "displacement_absolute_tolerance": 1.0e-9,
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
            "prebuild_dynamics"          : false,
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
            "body_domain_sub_model_part_list": [""]
        }""")

        this_defaults.AddMissingParameters(super().GetDefaultParameters())
        return this_defaults

    def PrepareModelPart(self):
        super().PrepareModelPart()
        KratosMultiphysics.Logger.PrintInfo("GeoMechanics_U_Pw_Solver", "Model reading finished.")

    def AddDofs(self):
        ## displacement dofs
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_X, KratosMultiphysics.REACTION_X,self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Y, KratosMultiphysics.REACTION_Y,self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Z, KratosMultiphysics.REACTION_Z,self.main_model_part)

        for node in self.main_model_part.Nodes:
            # adding TOTAL_DISPLACEMENT as dofs
            node.AddDof(KratosGeo.TOTAL_DISPLACEMENT_X)
            node.AddDof(KratosGeo.TOTAL_DISPLACEMENT_Y)
            node.AddDof(KratosGeo.TOTAL_DISPLACEMENT_Z)

        ## Fluid dofs
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.WATER_PRESSURE, KratosMultiphysics.REACTION_WATER_PRESSURE,self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosGeo.DT_WATER_PRESSURE, self.main_model_part)

        if self.settings["rotation_dofs"].GetBool():
            KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ROTATION_X, KratosMultiphysics.REACTION_MOMENT_X,self.main_model_part)
            KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ROTATION_Y, KratosMultiphysics.REACTION_MOMENT_Y,self.main_model_part)
            KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ROTATION_Z, KratosMultiphysics.REACTION_MOMENT_Z,self.main_model_part)

            KratosMultiphysics.VariableUtils().AddDof(KratosGeo.TOTAL_ROTATION_X, self.main_model_part)
            KratosMultiphysics.VariableUtils().AddDof(KratosGeo.TOTAL_ROTATION_Y, self.main_model_part)
            KratosMultiphysics.VariableUtils().AddDof(KratosGeo.TOTAL_ROTATION_Z, self.main_model_part)

        if (self.settings["solution_type"].GetString() == "Dynamic"):
            KratosMultiphysics.Logger.PrintInfo("GeoMechanics_U_Pw_Solver", "Dynamic analysis.")
            for node in self.main_model_part.Nodes:
                # adding VELOCITY as dofs
                node.AddDof(KratosMultiphysics.VELOCITY_X)
                node.AddDof(KratosMultiphysics.VELOCITY_Y)
                node.AddDof(KratosMultiphysics.VELOCITY_Z)
                # adding ACCELERATION as dofs
                node.AddDof(KratosMultiphysics.ACCELERATION_X)
                node.AddDof(KratosMultiphysics.ACCELERATION_Y)
                node.AddDof(KratosMultiphysics.ACCELERATION_Z)
                if self.settings["rotation_dofs"].GetBool():
                    # adding ANGULAR_VELOCITY as dofs
                    node.AddDof(KratosMultiphysics.ANGULAR_VELOCITY_X)
                    node.AddDof(KratosMultiphysics.ANGULAR_VELOCITY_Y)
                    node.AddDof(KratosMultiphysics.ANGULAR_VELOCITY_Z)
                    # adding ANGULAR_ACCELERATION as dofs
                    node.AddDof(KratosMultiphysics.ANGULAR_ACCELERATION_X)
                    node.AddDof(KratosMultiphysics.ANGULAR_ACCELERATION_Y)
                    node.AddDof(KratosMultiphysics.ANGULAR_ACCELERATION_Z)

        KratosMultiphysics.Logger.PrintInfo("GeoMechanics_U_Pw_Solver", "DOFs added correctly.")

    def Initialize(self):
        KratosMultiphysics.Logger.PrintInfo("::[GeoMechanics_U_Pw_Solver]:: ", "Initialisation ...")

        super().Initialize()

        KratosMultiphysics.Logger.PrintInfo("GeoMechanics_U_Pw_Solver", "solver.Initialize is set successfully")

        # Check if everything is assigned correctly
        self.Check()

        KratosMultiphysics.Logger.PrintInfo("GeoMechanics_U_Pw_Solver", "Solver initialization finished.")

    #### Specific internal functions ####

    def _ConstructScheme(self, scheme_type, solution_type):

        self.main_model_part.ProcessInfo.SetValue(KratosGeo.VELOCITY_COEFFICIENT,    1.0)
        self.main_model_part.ProcessInfo.SetValue(KratosGeo.DT_PRESSURE_COEFFICIENT, 1.0)

        if (scheme_type.lower() == "newmark"):
            beta  = self.settings["newmark_beta"].GetDouble()
            gamma = self.settings["newmark_gamma"].GetDouble()
            theta = self.settings["newmark_theta"].GetDouble()

            rayleigh_m = self.settings["rayleigh_m"].GetDouble()
            rayleigh_k = self.settings["rayleigh_k"].GetDouble()
            self.main_model_part.ProcessInfo.SetValue(KratosStructure.RAYLEIGH_ALPHA, rayleigh_m)
            self.main_model_part.ProcessInfo.SetValue(KratosStructure.RAYLEIGH_BETA,  rayleigh_k)

            KratosMultiphysics.Logger.PrintInfo("GeoMechanics_U_Pw_Solver, solution_type", solution_type)
            if (solution_type.lower() == "quasi-static" or solution_type.lower() == "quasi_static"):
                if (rayleigh_m < 1.0e-20 and rayleigh_k < 1.0e-20):
                    KratosMultiphysics.Logger.PrintInfo("GeoMechanics_U_Pw_Solver, scheme", "Quasi-UnDamped.")
                    return KratosGeo.NewmarkQuasistaticUPwScheme(beta, gamma, theta)
                else:
                    KratosMultiphysics.Logger.PrintInfo("GeoMechanics_U_Pw_Solver, scheme", "Quasi-Damped.")
                    return KratosGeo.NewmarkQuasistaticDampedUPwScheme(beta, gamma, theta)

            if solution_type.lower() == "dynamic":
                KratosMultiphysics.Logger.PrintInfo("GeoMechanics_U_Pw_Solver, scheme", "Dynamic.")
                return KratosGeo.NewmarkDynamicUPwScheme(beta, gamma, theta)

            raise RuntimeError(f"Undefined solution type '{solution_type}'")

        if (scheme_type.lower() == "backward_euler" or scheme_type.lower() == "backward-euler"):
            if (solution_type.lower() == "quasi-static" or solution_type.lower() == "quasi_static"):
                KratosMultiphysics.Logger.PrintInfo("GeoMechanics_U_Pw_Solver, scheme", "Backward Euler.")
                return KratosGeo.BackwardEulerQuasistaticUPwScheme()

            raise RuntimeError(f"Undefined/incompatible solution type with Backward Euler: '{solution_type}'")

        raise RuntimeError("Apart from Newmark and Backward Euler, other scheme_type are not available.")

    def _ConstructConvergenceCriterion(self, convergence_criterion):

        if convergence_criterion.lower() == "displacement_criterion":
            return self._MakeDisplacementCriterion()

        if convergence_criterion.lower() == "residual_criterion":
           return self._MakeResidualCriterion()

        if convergence_criterion.lower() == "and_criterion":
            return KratosMultiphysics.AndCriteria(self._MakeResidualCriterion(), self._MakeDisplacementCriterion())

        if convergence_criterion.lower() == "or_criterion":
            return KratosMultiphysics.OrCriteria(self._MakeResidualCriterion(),  self._MakeDisplacementCriterion())

        if convergence_criterion.lower() == "water_pressure_criterion":
            return self._MakeWaterPressureCriterion()

        if convergence_criterion.lower() == "displacement_and_water_pressure_criterion":
            d_and_pw_convergence_criterion = KratosMultiphysics.MixedGenericCriteria([(self._MakeDisplacementCriterion()),(self._MakeWaterPressureCriterion())])
            d_and_pw_convergence_criterion.SetEchoLevel(self.settings["echo_level"].GetInt())
            return d_and_pw_convergence_criterion

        err_msg =  "The requested convergence criterion \"" + convergence_criterion + "\" is not available!\n"
        err_msg += "Available options are: \"displacement_criterion\", \"residual_criterion\", \"and_criterion\", \"or_criterion\", \"water_pressure_criterion\", \"displacement_and_water_pressure_criterion\""
        raise RuntimeError(err_msg)

    def _MakeDisplacementCriterion(self):
        relative_tolerance = self.settings["displacement_relative_tolerance"].GetDouble()
        absolute_tolerance = self.settings["displacement_absolute_tolerance"].GetDouble()

        displacement_criterion = KratosMultiphysics.DisplacementCriteria(relative_tolerance, absolute_tolerance)
        displacement_criterion.SetEchoLevel(self.settings["echo_level"].GetInt())

        return displacement_criterion

    def _CreateBuilderAndSolver(self):
        if (self.settings["block_builder"].GetBool() and
            self.settings.Has("prebuild_dynamics") and self.settings["prebuild_dynamics"].GetBool()):
            return KratosGeo.ResidualBasedBlockBuilderAndSolverWithMassAndDamping(self.linear_solver)

        return super()._CreateBuilderAndSolver()

    def _CheckConvergence(self):
        return self.solving_strategy.IsConverged()
