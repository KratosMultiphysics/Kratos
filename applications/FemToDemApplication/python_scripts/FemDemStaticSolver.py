from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.SolidMechanicsApplication as KratosSolid
import KratosMultiphysics.FemToDemApplication as KratosFemDem

# Import the mechanical solver base class
import KratosMultiphysics.FemToDemApplication.FemDemMechanicalSolver as BaseSolver

def CreateSolver(main_model_part, custom_settings):
    return StaticMechanicalSolver(main_model_part, custom_settings)

class StaticMechanicalSolver(BaseSolver.FemDemMechanicalSolver):

    """The solid mechanics static solver.
    This class creates the mechanical solvers for static analysis.
    Public member variables:
    See solid_mechanics_solver.py for more information.
    """
    def __init__(self, main_model_part, custom_settings):

        # Set defaults and validate custom settings.
        static_settings = KratosMultiphysics.Parameters("""
        {
        }
        """)
        self._validate_and_transfer_matching_settings(custom_settings, static_settings)
        # Validate the remaining settings in the base class.
        if not custom_settings.Has("solution_type"):
            custom_settings.AddEmptyValue("solution_type")
            custom_settings["solution_type"].SetString("Static") # Override defaults in the base class.
        if not custom_settings.Has("scheme_type"):
            custom_settings.AddEmptyValue("scheme_type")
            custom_settings["scheme_type"].SetString("Non-Linear") # Override defaults in the base class.
        if not custom_settings.Has("extrapolation_required"):
            custom_settings.AddEmptyValue("extrapolation_required")
            custom_settings["extrapolation_required"].SetBool(False)

        # Construct the base solver.
        super(StaticMechanicalSolver, self).__init__(main_model_part, custom_settings)

        print("::[Static_Mechanical_Solver]:: Constructed")

    #### Solver internal methods ####

    def _create_solution_scheme (self):
        scheme_type = self.settings["scheme_type"].GetString()
        if(scheme_type == "Linear"):
            mechanical_scheme = KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme()
        elif(scheme_type == "Non-Linear" ):
            if(self.settings["component_wise"].GetBool() == True):
                dynamic_factor = 0.0
                damp_factor_m  = 0.0
                mechanical_scheme = KratosSolid.ComponentWiseBossakScheme(damp_factor_m, dynamic_factor)
            else:
                mechanical_scheme = KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme()
        elif(scheme_type == "RotationStatic"):
            dynamic_factor = 0.0
            damp_factor_m  = 0.0
            mechanical_scheme = KratosSolid.ResidualBasedRotationNewmarkScheme(dynamic_factor, damp_factor_m)
        elif(scheme_type == "RotationEMC"):
            dynamic_factor = 0.0
            mechanical_scheme = KratosSolid.ResidualBasedRotationEMCScheme(dynamic_factor)
        else:
            raise Exception("Unsupported scheme_type: " + scheme_type)

        return mechanical_scheme

    def _create_mechanical_solver(self):
        if(self.settings["component_wise"].GetBool() == True):
            mechanical_solver = self._create_component_wise_strategy()
        elif(self.settings["line_search"].GetBool() == True):
            if(self.settings["implex"].GetBool() == True):
                mechanical_solver = self._create_line_search_implex_strategy()
            else:
                mechanical_solver = self._create_line_search_strategy()
        else:
            if(self.settings["scheme_type"].GetString() == "Linear"):
                mechanical_solver = self._create_linear_strategy()
            else:
                if self.settings["extrapolation_required"].GetBool():
                    mechanical_solver = self._create_newton_raphson_hexaedrons_strategy()
                else:
                    mechanical_solver = self._create_newton_raphson_strategy()

        return mechanical_solver


    def _create_component_wise_strategy(self):
        computing_model_part = self.GetComputingModelPart()
        mechanical_scheme = self._get_solution_scheme()
        linear_solver = self._get_linear_solver()
        mechanical_convergence_criterion = self._get_convergence_criterion()
        builder_and_solver = self._get_builder_and_solver()
        return KratosSolid.ComponentWiseNewtonRaphsonStrategy(computing_model_part,
                                                              mechanical_scheme,
                                                              linear_solver,
                                                              mechanical_convergence_criterion,
                                                              builder_and_solver,
                                                              self.settings["max_iteration"].GetInt(),
                                                              self.settings["compute_reactions"].GetBool(),
                                                              self.settings["reform_dofs_at_each_step"].GetBool(),
                                                              self.settings["move_mesh_flag"].GetBool())

    def _create_line_search_strategy(self):
        computing_model_part = self.GetComputingModelPart()
        mechanical_scheme = self._get_solution_scheme()
        linear_solver = self._get_linear_solver()
        mechanical_convergence_criterion = self._get_convergence_criterion()
        builder_and_solver = self._get_builder_and_solver()
        return KratosSolid.ResidualBasedNewtonRaphsonLineSearchStrategy(computing_model_part,
                                                                        mechanical_scheme,
                                                                        linear_solver,
                                                                        mechanical_convergence_criterion,
                                                                        builder_and_solver,
                                                                        self.settings["max_iteration"].GetInt(),
                                                                        self.settings["compute_reactions"].GetBool(),
                                                                        self.settings["reform_dofs_at_each_step"].GetBool(),
                                                                        self.settings["move_mesh_flag"].GetBool())

    def _create_line_search_implex_strategy(self):
        computing_model_part = self.GetComputingModelPart()
        mechanical_scheme = self._get_solution_scheme()
        linear_solver = self._get_linear_solver()
        mechanical_convergence_criterion = self._get_convergence_criterion()
        builder_and_solver = self._get_builder_and_solver()
        return KratosSolid.ResidualBasedNewtonRaphsonLineSearchImplexStrategy(computing_model_part,
                                                                              mechanical_scheme,
                                                                              linear_solver,
                                                                              mechanical_convergence_criterion,
                                                                              builder_and_solver,
                                                                              self.settings["max_iteration"].GetInt(),
                                                                              self.settings["compute_reactions"].GetBool(),
                                                                              self.settings["reform_dofs_at_each_step"].GetBool(),
                                                                              self.settings["move_mesh_flag"].GetBool())


    def _create_linear_strategy(self):
        computing_model_part = self.GetComputingModelPart()
        mechanical_scheme = self._get_solution_scheme()
        linear_solver = self._get_linear_solver()
        builder_and_solver = self._get_builder_and_solver()
        return KratosMultiphysics.ResidualBasedLinearStrategy(computing_model_part,
                                                              mechanical_scheme,
                                                              linear_solver,
                                                              builder_and_solver,
                                                              self.settings["compute_reactions"].GetBool(),
                                                              self.settings["reform_dofs_at_each_step"].GetBool(),
                                                              False,
                                                              self.settings["move_mesh_flag"].GetBool())


    def _create_newton_raphson_strategy(self):
        computing_model_part = self.GetComputingModelPart()
        mechanical_scheme = self._get_solution_scheme()
        linear_solver = self._get_linear_solver()
        mechanical_convergence_criterion = self._get_convergence_criterion()
        builder_and_solver = self._get_builder_and_solver()
        return KratosMultiphysics.ResidualBasedNewtonRaphsonStrategy(computing_model_part,
                                                                     mechanical_scheme,
                                                                     linear_solver,
                                                                     mechanical_convergence_criterion,
                                                                     builder_and_solver,
                                                                     self.settings["max_iteration"].GetInt(),
                                                                     self.settings["compute_reactions"].GetBool(),
                                                                     self.settings["reform_dofs_at_each_step"].GetBool(),
                                                                     self.settings["move_mesh_flag"].GetBool())

    def _create_newton_raphson_hexaedrons_strategy(self):
        computing_model_part = self.GetComputingModelPart()
        mechanical_scheme = self._get_solution_scheme()
        linear_solver = self._get_linear_solver()
        mechanical_convergence_criterion = self._get_convergence_criterion()
        builder_and_solver = self._get_builder_and_solver()
        return KratosFemDem.HexahedraNewtonRaphsonStrategy(computing_model_part,
                                                                     mechanical_scheme,
                                                                     linear_solver,
                                                                     mechanical_convergence_criterion,
                                                                     self.settings["max_iteration"].GetInt(),
                                                                     self.settings["compute_reactions"].GetBool(),
                                                                     self.settings["reform_dofs_at_each_step"].GetBool(),
                                                                     self.settings["move_mesh_flag"].GetBool())

