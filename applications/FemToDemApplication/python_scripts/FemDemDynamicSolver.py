#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.FemToDemApplication as KratosFemDem
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication

# Import the mechanical solver base class
import KratosMultiphysics.FemToDemApplication.FemDemMechanicalSolver as BaseSolver

def KratosPrintInfo(message):
    KratosMultiphysics.Logger.Print(message, label="")
    KratosMultiphysics.Logger.Flush()

def CreateSolver(main_model_part, custom_settings, DEMStrategy = None):
    return ImplicitMechanicalSolver(main_model_part, custom_settings, DEMStrategy)

class ImplicitMechanicalSolver(BaseSolver.FemDemMechanicalSolver):

    def __init__(self, main_model_part, custom_settings, DEMStrategy = None):

        # Set defaults and validate custom settings.
        self.dynamic_settings = KratosMultiphysics.Parameters("""
        {
            "damp_factor_m" :-0.01,
            "dynamic_factor": 1.0,
            "rayleigh_damping": false,
            "rayleigh_alpha": 0.0,
            "rayleigh_beta" : 0.0
        }
        """)

        self._validate_and_transfer_matching_settings(custom_settings, self.dynamic_settings)
        # Validate the remaining settings in the base class.
        if not custom_settings.Has("scheme_type"): # Override defaults in the base class.
            custom_settings.AddEmptyValue("scheme_type")
            custom_settings["scheme_type"].SetString("Newmark")

        if not custom_settings.Has("extrapolation_required"):
            custom_settings.AddEmptyValue("extrapolation_required")
            custom_settings["extrapolation_required"].SetBool(False)

        # Construct the base solver.
        super(ImplicitMechanicalSolver, self).__init__(main_model_part, custom_settings)

        print("::[Implicit_Dynamic_Solver]:: Constructed")

        self.DEMStrategy = None
        if DEMStrategy != None:
            self.DEMStrategy = DEMStrategy

    #### Solver internal methods ####

    def _create_solution_scheme(self):

        scheme_type = self.settings["scheme_type"].GetString()

        if( self.dynamic_settings["rayleigh_damping"].GetBool() == True ):
            self.main_model_part.ProcessInfo[KratosFemDem.RAYLEIGH_ALPHA] = self.dynamic_settings["rayleigh_alpha"].GetDouble()
            self.main_model_part.ProcessInfo[KratosFemDem.RAYLEIGH_BETA]  = self.dynamic_settings["rayleigh_beta"].GetDouble()
        else:
            self.main_model_part.ProcessInfo[KratosFemDem.RAYLEIGH_ALPHA] = 0.0
            self.main_model_part.ProcessInfo[KratosFemDem.RAYLEIGH_BETA]  = 0.0

        if(scheme_type == "Newmark"):
            damp_factor_m = 0.0
            mechanical_scheme = KratosMultiphysics.ResidualBasedBossakDisplacementScheme(damp_factor_m)
        elif(scheme_type == "Bossak"):
            damp_factor_m = self.dynamic_settings["damp_factor_m"].GetDouble()
            mechanical_scheme = KratosMultiphysics.ResidualBasedBossakDisplacementScheme(damp_factor_m)
        elif(scheme_type == "central_differences"):
            mechanical_scheme = StructuralMechanicsApplication.ExplicitCentralDifferencesScheme(self.settings["max_delta_time"].GetDouble(),
                                                                             self.settings["fraction_delta_time"].GetDouble(),
                                                                             self.settings["time_step_prediction_level"].GetDouble())
        else:
            raise Exception("Unsupported scheme_type: " + scheme_type)

        return mechanical_scheme

    def _create_mechanical_solver(self):
        if (self.settings["line_search"].GetBool()):
            mechanical_solver = self._create_line_search_strategy()
        else:
            if self.settings["extrapolation_required"].GetBool():
                mechanical_solver = self._create_newton_raphson_hexaedrons_strategy()
            elif (self.settings["strategy_type"] == "arc_length"):
                mechanical_solver = self._create_ramm_arc_length_strategy()
            else:
                if self.DEMStrategy == None:
                    if self.settings["time_integration_method"].GetString() == "Implicit":
                        mechanical_solver = self._create_newton_raphson_strategy()
                    else:
                        KratosPrintInfo("  ______            _ _      _ _   " + "\n"
                                        " |  ____|          | (_)    (_) |  " +"\n"
                                        " | |__  __  ___ __ | |_  ___ _| |_ " +"\n"
                                        " |  __| \ \/ / '_ \| | |/ __| | __|" +"\n"
                                        " | |____ >  <| |_) | | | (__| | |_ " +"\n"
                                        " |______/_/\_\ .__/|_|_|\___|_|\__| FEM-DEM version" +"\n"
                                        "             | |                   " +"\n"
                                        "             |_|                   ")
                        mechanical_solver = self._create_explicit_strategy()
                else:
                    mechanical_solver = self._create_DEM_coupled_newton_raphson_strategy(self.DEMStrategy)
        return mechanical_solver

    def _create_line_search_strategy(self):
        computing_model_part = self.GetComputingModelPart()
        mechanical_scheme = self._get_solution_scheme()
        linear_solver = self._get_linear_solver()
        mechanical_convergence_criterion = self._get_convergence_criterion()
        builder_and_solver = self._get_builder_and_solver()
        return KratosMultiphysics.LineSearchStrategy(computing_model_part,
                                                     mechanical_scheme,
                                                     linear_solver,
                                                     mechanical_convergence_criterion,
                                                     builder_and_solver,
                                                     self.settings["max_iteration"].GetInt(),
                                                     self.settings["compute_reactions"].GetBool(),
                                                     self.settings["reform_dofs_at_each_step"].GetBool(),
                                                     self.settings["move_mesh_flag"].GetBool())

    def _create_newton_raphson_strategy(self):
        computing_model_part = self.GetComputingModelPart()
        mechanical_scheme = self._get_solution_scheme()
        mechanical_convergence_criterion = self._get_convergence_criterion()
        builder_and_solver = self._get_builder_and_solver()
        return KratosMultiphysics.ResidualBasedNewtonRaphsonStrategy(computing_model_part,
                                                                     mechanical_scheme,
                                                                     mechanical_convergence_criterion,
                                                                     builder_and_solver,
                                                                     self.settings["max_iteration"].GetInt(),
                                                                     self.settings["compute_reactions"].GetBool(),
                                                                     self.settings["reform_dofs_at_each_step"].GetBool(),
                                                                     self.settings["move_mesh_flag"].GetBool())

    def _create_explicit_strategy(self):
        computing_model_part = self.GetComputingModelPart()
        mechanical_scheme = self._get_solution_scheme()

        return StructuralMechanicsApplication.MechanicalExplicitStrategy(computing_model_part,
                                            mechanical_scheme,
                                            self.settings["compute_reactions"].GetBool(),
                                            self.settings["reform_dofs_at_each_step"].GetBool(),
                                            self.settings["move_mesh_flag"].GetBool())

    def _create_DEM_coupled_newton_raphson_strategy(self, DEM_strategy):
        computing_model_part = self.GetComputingModelPart()
        mechanical_scheme = self._get_solution_scheme()
        mechanical_convergence_criterion = self._get_convergence_criterion()
        builder_and_solver = self._get_builder_and_solver()
        return KratosFemDem.ResidualBasedDEMCoupledNewtonRaphsonStrategy(computing_model_part,
                                                                         DEM_strategy.cplusplus_strategy,
                                                                         mechanical_scheme,
                                                                         mechanical_convergence_criterion,
                                                                         builder_and_solver,
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
        for i in range(self.loads_sub_sub_model_part_list.size()):
            computing_sub_model = computing_model_part.CreateSubModelPart("sub_" + self.settings["arc_length_loads_sub_model_part_list"][i].GetString())
            for node in self.main_model_part.GetSubModelPart(self.settings["arc_length_loads_sub_model_part_list"][i].GetString()).Nodes:
                computing_sub_model.AddNode(node, 0)

        computing_model_part = self.GetComputingModelPart()
        mechanical_scheme = self._get_solution_scheme()
        linear_solver = self._get_linear_solver()
        mechanical_convergence_criterion = self._get_convergence_criterion()
        builder_and_solver = self._get_builder_and_solver()

        # Arc-Length strategy
        self.main_model_part.ProcessInfo.SetValue(KratosPoro.ARC_LENGTH_LAMBDA,1.0)
        self.main_model_part.ProcessInfo.SetValue(KratosPoro.ARC_LENGTH_RADIUS_FACTOR,1.0)

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