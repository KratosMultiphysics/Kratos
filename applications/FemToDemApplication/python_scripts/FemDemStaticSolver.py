import KratosMultiphysics
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
        if (self.settings["component_wise"].GetBool()):
            mechanical_solver = self._create_component_wise_strategy()
        elif (self.settings["line_search"].GetBool()):
            if(self.settings["implex"].GetBool()):
                mechanical_solver = self._create_line_search_implex_strategy()
            else:
                mechanical_solver = self._create_line_search_strategy()
        else:
            if(self.settings["scheme_type"].GetString() == "Linear"):
                mechanical_solver = self._create_linear_strategy()
            else:
                if self.settings["extrapolation_required"].GetBool():
                    mechanical_solver = self._create_newton_raphson_hexaedrons_strategy()
                elif (self.settings["strategy_type"].GetString() == "arc_length"):
                    mechanical_solver = self._create_ramm_arc_length_strategy()
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

    def _create_ramm_arc_length_strategy(self):
        # Create list of sub sub model parts (it is a copy of the standard lists with a different name)
        import json
        self.loads_sub_sub_model_part_list = []
        for i in range(self.settings["arc_length_loads_sub_model_part_list"].size()):
            self.loads_sub_sub_model_part_list.append("sub_" + self.settings["arc_length_loads_sub_model_part_list"][i].GetString())
        self.loads_sub_sub_model_part_list = KratosMultiphysics.Parameters(json.dumps(self.loads_sub_sub_model_part_list))

        computing_model_part = self.GetComputingModelPart()

        # We create the arc-length loads subsubmodel parts
        # Adding the nodes and conditions
        for i in range(self.loads_sub_sub_model_part_list.size()):
            computing_sub_model = computing_model_part.CreateSubModelPart("sub_" + self.settings["arc_length_loads_sub_model_part_list"][i].GetString())
            sub_model_part = self.main_model_part.GetSubModelPart(self.settings["arc_length_loads_sub_model_part_list"][i].GetString())
            for node in sub_model_part.Nodes:
                computing_sub_model.AddNode(node, 0)
            for cond in sub_model_part.Conditions:
                computing_sub_model.AddCondition(cond, 0)

        mechanical_scheme = self._get_solution_scheme()
        linear_solver = self._get_linear_solver()
        mechanical_convergence_criterion = self._get_convergence_criterion()
        builder_and_solver = self._get_builder_and_solver()

        # Arc-Length strategy
        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.ARC_LENGTH_LAMBDA,1.0)
        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.ARC_LENGTH_RADIUS_FACTOR,1.0)

        self.strategy_params = KratosMultiphysics.Parameters("{}")
        self.strategy_params.AddValue("desired_iterations", self.settings["arc_length_desired_iterations"])
        self.strategy_params.AddValue("max_radius_factor", self.settings["arc_length_max_radius_factor"])
        self.strategy_params.AddValue("min_radius_factor", self.settings["arc_length_min_radius_factor"])
        self.strategy_params.AddValue("loads_sub_model_part_list", self.loads_sub_sub_model_part_list)
        self.strategy_params.AddValue("loads_variable_list", self.settings["arc_length_loads_variable_list"])

        return KratosMultiphysics.RammArcLengthStrategy(computing_model_part,
                                                        mechanical_scheme,
                                                        linear_solver,
                                                        mechanical_convergence_criterion,
                                                        builder_and_solver,
                                                        self.strategy_params,
                                                        self.settings["max_iteration"].GetInt(),
                                                        self.settings["compute_reactions"].GetBool(),
                                                        self.settings["reform_dofs_at_each_step"].GetBool(),
                                                        self.settings["move_mesh_flag"].GetBool())