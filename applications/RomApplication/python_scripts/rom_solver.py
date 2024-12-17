# Importing the Kratos Library
from typing import List
import KratosMultiphysics

# Import applications
import KratosMultiphysics.RomApplication as KratosROM

def CreateSolver(cls, model, custom_settings):
    class ROMSolver(cls):
        """ROM solver generic main class.
        This class serves as a generic class to make a standard Kratos solver ROM-compatible.
        It extends the default parameters to include the \'rom_settings\' and overrides the
        creation of the builder and solver to use the ROM one.
        """

        def __init__(self, model, custom_settings):
            super().__init__(model, custom_settings)
            KratosMultiphysics.Logger.PrintInfo("::[ROMSolver]:: ", "Construction finished")

            this_analysis_stage = self.settings["rom_settings"]["analysis_stage"].GetString() #TODO Could we simply use the model part name or other settings as proxi here?

            # if this_analysis_stage == "KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis":
            #     self._CreateSolutionStrategy = self.Rom_Structural_CreateSolutionStrategy
            # elif this_analysis_stage == "KratosMultiphysics.FluidDynamicsApplication.fluid_dynamics_analysis":
            #     self._CreateSolutionStrategy = self.Rom_Fluid_CreateSolutionStrategy
            # elif this_analysis_stage == "KratosMultiphysics.CompressiblePotentialFlowApplication.potential_flow_analysis":
            #     self._CreateSolutionStrategy = self.Rom_CompressiblePotential_CreateSolutionStrategy
            # else:
            #     error
            # # elif this_analysis_stage == "KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis":
            # #     pass
            # # elif this_analysis_stage == "KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis":
            # #     pass

            # #eliminate from parameters so that it doesnt botther the B&S
            self.settings["rom_settings"].RemoveValue("analysis_stage")



        @classmethod
        def GetDefaultParameters(cls):
            default_settings = KratosMultiphysics.Parameters("""{
                "projection_strategy" : "galerkin",
                "assembling_strategy" : "global",
                "rom_settings": {
                    "nodal_unknowns": [],
                    "number_of_rom_dofs": 0,
                    "rom_bns_settings": {}
                }
            }""")
            default_settings.AddMissingParameters(super().GetDefaultParameters())
            return default_settings

        def _CreateBuilderAndSolver(self):
            linear_solver = self._GetLinearSolver()
            rom_parameters, solving_strategy = self._ValidateAndReturnRomParameters()
            available_solving_strategies = {
                "elemental_galerkin": KratosROM.ROMBuilderAndSolver,
                "global_galerkin": KratosROM.GlobalROMBuilderAndSolver,
                "lspg": KratosROM.LeastSquaresPetrovGalerkinROMBuilderAndSolver,
                "lspg_ann": KratosROM.AnnPromLeastSquaresPetrovGalerkinROMBuilderAndSolver,
                "galerkin_ann": KratosROM.AnnPromGlobalROMBuilderAndSolver,
                "elemental_petrov_galerkin": KratosROM.PetrovGalerkinROMBuilderAndSolver,
                "global_petrov_galerkin": KratosROM.GlobalPetrovGalerkinROMBuilderAndSolver
            }
            if solving_strategy in available_solving_strategies:
                return available_solving_strategies[solving_strategy](linear_solver, rom_parameters)
            else:
                err_msg = f"'Solving_strategy': '{solving_strategy}' is not available. Please select one of the following: {list(available_solving_strategies.keys())}."
                raise ValueError(err_msg)

        def _ValidateAndReturnRomParameters(self):
            # Check that the number of ROM DOFs has been provided
            n_rom_dofs = self.settings["rom_settings"]["number_of_rom_dofs"].GetInt()
            if not n_rom_dofs > 0:
                err_msg = "\'number_of_rom_dofs\' in \'rom_settings\' is {}. Please set a larger than zero value.".format(n_rom_dofs)
                raise Exception(err_msg)

            # Check if the nodal unknowns have been provided by the user
            # If not, take the DOFs list from the base solver
            nodal_unknowns = self.settings["rom_settings"]["nodal_unknowns"].GetStringArray()
            if len(nodal_unknowns) == 0:
                solver_dofs_list = self.GetDofsList()
                if not len(solver_dofs_list) == 0:
                    self.settings["rom_settings"]["nodal_unknowns"].SetStringArray(solver_dofs_list)
                else:
                    err_msg = "\'nodal_unknowns\' in \'rom_settings\' is not provided and there is a not-valid implementation in base solver."
                    err_msg += " Please manually set \'nodal_unknowns\' in \'rom_settings\'."
                    raise Exception(err_msg)
            projection_strategy = self.settings["projection_strategy"].GetString()
            assembling_strategy = self.settings["assembling_strategy"].GetString()

            # Check for the 'elemental' assembling strategy in case of 'lspg' projection strategy
            if projection_strategy == "lspg" and assembling_strategy == "elemental":
                warn_msg = "Elemental LSPG is not yet available. Using default global assembling strategy instead."
                KratosMultiphysics.Logger.PrintWarning("::[ROMSolver]:: ", warn_msg)

            # For now, only Galerkin and Petrov-Galerkin projections have the elemental or global approach option
            if projection_strategy in ("galerkin", "petrov_galerkin"): #TODO: Possibility of doing elemental lspg
                available_assembling_strategies = {
                    "global",
                    "elemental"
                }
                if assembling_strategy in available_assembling_strategies:
                    # Add the assembling strategy prefix to the projection strategy
                    projection_strategy = f"{assembling_strategy}_{projection_strategy}"
                else:
                    err_msg = f"'Assembling_strategy': '{assembling_strategy}' is not available. Please select one of the following: {list(available_assembling_strategies)}."
                    raise ValueError(err_msg)

            self._AssignMissingInnerRomParameters(projection_strategy)

            # Check that the number of ROM DOFs has been provided
            if projection_strategy in ("galerkin_ann", "lspg_ann") == False:
                n_rom_dofs = self.settings["rom_settings"]["number_of_rom_dofs"].GetInt()
                if not n_rom_dofs > 0:
                    err_msg = "\'number_of_rom_dofs\' in \'rom_settings\' is {}. Please set a larger than zero value.".format(n_rom_dofs)
                    raise Exception(err_msg)                

            # Return the validated ROM parameters
            return self.settings["rom_settings"], projection_strategy

        def _AssignMissingInnerRomParameters(self, projection_strategy):
            if not self.settings["rom_settings"].Has("rom_bns_settings"):
                self.settings["rom_settings"].AddEmptyValue("rom_bns_settings")
            monotonicity_preserving = self.settings["rom_settings"]["rom_bns_settings"]["monotonicity_preserving"].GetBool() if self.settings["rom_settings"]["rom_bns_settings"].Has("monotonicity_preserving") else False
            if projection_strategy in ("global_galerkin", "lspg", "global_petrov_galerkin", "galerkin_ann", "lspg_ann"):
                self.settings["rom_settings"]["rom_bns_settings"].AddBool("monotonicity_preserving", monotonicity_preserving)

        # #####################################################
        # #  Setting strategy for Structural Mechanics App   #
        # #####################################################

        # def Rom_Structural_CreateSolutionStrategy(self):
        #     analysis_type = self.settings["analysis_type"].GetString()
        #     if analysis_type == "linear":
        #         mechanical_solution_strategy = self._create_linear_strategy() #TODO should this stay as is
        #     elif analysis_type == "non_linear":
        #         # Create strategy
        #         if self.settings["solving_strategy_settings"]["type"].GetString() == "newton_raphson":
        #             mechanical_solution_strategy = self._structural_create_newton_raphson_strategy_ROM()
        #         elif self.settings["solving_strategy_settings"]["type"].GetString() == "line_search":
        #             error #TODO implement
        #             #mechanical_solution_strategy = self._create_line_search_strategy_ROM()
        #         elif self.settings["solving_strategy_settings"]["type"].GetString() == "arc_length":
        #             error #TODO implement
        #             #mechanical_solution_strategy = self._create_arc_length_strategy()
        #     else:
        #         err_msg =  "The requested analysis type \"" + analysis_type + "\" is not available!\n"
        #         err_msg += "Available options are: \"linear\", \"non_linear\""
        #         raise Exception(err_msg)
        #     return mechanical_solution_strategy


        # def _structural_create_newton_raphson_strategy_ROM(self):
        #     computing_model_part = self.GetComputingModelPart()
        #     mechanical_scheme = self._GetScheme()
        #     mechanical_convergence_criterion = self._GetConvergenceCriterion()
        #     builder_and_solver = self._GetBuilderAndSolver()
        #     #this calls the ROM strategy
        #     strategy = KratosMultiphysics.RomApplication.RomResidualBasedNewtonRaphsonStrategy(computing_model_part,
        #                                                                 mechanical_scheme,
        #                                                                 mechanical_convergence_criterion,
        #                                                                 builder_and_solver,
        #                                                                 self.settings["max_iteration"].GetInt(),
        #                                                                 self.settings["compute_reactions"].GetBool(),
        #                                                                 self.settings["reform_dofs_at_each_step"].GetBool(),
        #                                                                 self.settings["move_mesh_flag"].GetBool())
        #     strategy.SetUseOldStiffnessInFirstIterationFlag(self.settings["use_old_stiffness_in_first_iteration"].GetBool())
        #     return strategy



        # #####################################################
        # #      Setting strategy for Fluid Dynamics App      #
        # #####################################################

        # def Rom_Fluid_CreateSolutionStrategy(self):
        #     analysis_type = self.settings["analysis_type"].GetString()
        #     if analysis_type == "linear":
        #         solution_strategy = self._CreateLinearStrategy() #TODO should this stay as is
        #     elif analysis_type == "non_linear":
        #         solution_strategy = self._fluid_CreateNewtonRaphsonStrategy_ROM()
        #     else:
        #         err_msg =  "The requested analysis type \"" + analysis_type + "\" is not available!\n"
        #         err_msg += "Available options are: \"linear\", \"non_linear\""
        #         raise Exception(err_msg)
        #     solution_strategy.SetEchoLevel(self.settings["echo_level"].GetInt())
        #     return solution_strategy


        # def _CreateLinearStrategy(self):
        #     computing_model_part = self.GetComputingModelPart()
        #     time_scheme = self._GetScheme()
        #     builder_and_solver = self._GetBuilderAndSolver()
        #     calculate_norm_dx = False
        #     return KratosMultiphysics.ResidualBasedLinearStrategy(
        #         computing_model_part,
        #         time_scheme,
        #         builder_and_solver,
        #         self.settings["compute_reactions"].GetBool(),
        #         self.settings["reform_dofs_at_each_step"].GetBool(),
        #         calculate_norm_dx,
        #         self.settings["move_mesh_flag"].GetBool())
    
        # def _fluid_CreateNewtonRaphsonStrategy_ROM(self):
        #     computing_model_part = self.GetComputingModelPart()
        #     time_scheme = self._GetScheme()
        #     convergence_criterion = self._GetConvergenceCriterion()
        #     builder_and_solver = self._GetBuilderAndSolver()
        #     return KratosMultiphysics.RomApplication.RomResidualBasedNewtonRaphsonStrategy(
        #         computing_model_part,
        #         time_scheme,
        #         convergence_criterion,
        #         builder_and_solver,
        #         self.settings["maximum_iterations"].GetInt(),
        #         self.settings["compute_reactions"].GetBool(),
        #         self.settings["reform_dofs_at_each_step"].GetBool(),
        #         self.settings["move_mesh_flag"].GetBool())
        

        # def _CreateLineSearchStrategy(self):
        #     if self.settings["solving_strategy_settings"].Has("advanced_settings"):
        #         settings = self.settings["solving_strategy_settings"]["advanced_settings"]
        #         settings.AddMissingParameters(self._GetDefaultLineSearchParameters())
        #     else:
        #         settings = self._GetDefaultLineSearchParameters()
        #     settings.AddValue("max_iteration", self.settings["maximum_iterations"])
        #     settings.AddValue("compute_reactions", self.settings["compute_reactions"])
        #     settings.AddValue("reform_dofs_at_each_step", self.settings["reform_dofs_at_each_step"])
        #     settings.AddValue("move_mesh_flag", self.settings["move_mesh_flag"])
        #     computing_model_part = self.GetComputingModelPart()
        #     time_scheme = self._GetScheme()
        #     convergence_criterion = self._GetConvergenceCriterion()
        #     builder_and_solver = self._GetBuilderAndSolver()
        #     solution_strategy = KratosROM.RomLineSearchStrategy(computing_model_part,
        #         time_scheme,
        #         convergence_criterion,
        #         builder_and_solver,
        #         settings)
        #     return solution_strategy
        
        # @classmethod
        # def _GetDefaultLineSearchParameters(self):
        #     default_line_search_parameters = KratosMultiphysics.Parameters(r"""{
        #             "max_line_search_iterations" : 5,
        #             "first_alpha_value"          : 0.5,
        #             "second_alpha_value"         : 1.0,
        #             "min_alpha"                  : 0.1,
        #             "max_alpha"                  : 1.0,
        #             "line_search_tolerance"      : 0.5
        #         }""")
        #     return default_line_search_parameters



        # ############################################################
        # #      Setting strategy for  Compressible Potential App    #
        # ############################################################

        # def _GetAnalysisType(self):
        #     if self.settings["formulation"].Has("element_type"):
        #         element_type = self.settings["formulation"]["element_type"].GetString()
        #         if (element_type == "perturbation_incompressible"):
        #             return "linear"
        #         elif (element_type == "perturbation_compressible"
        #             or element_type == "perturbation_transonic"):
        #             return "non_linear"
        #         else:
        #             raise RuntimeError(f"Argument {element_type} not found in formulation settings.")


        # def Rom_CompressiblePotential_CreateSolutionStrategy(self):
        #     analysis_type = self._GetAnalysisType()
        #     if analysis_type == "linear":
        #         # Create strategy
        #         solution_strategy = self._CreateLinearStrategy()
        #     elif analysis_type == "non_linear":
        #         # Create strategy
        #         if self.settings["solving_strategy_settings"]["type"].GetString() == "newton_raphson":
        #             solution_strategy = self._fluid_CreateNewtonRaphsonStrategy_ROM() #same as in regular fluid
        #         elif self.settings["solving_strategy_settings"]["type"].GetString() == "line_search":
        #             solution_strategy = self._CreateLineSearchStrategy()
        #     else:
        #         err_msg = f"Unknown strategy type: {analysis_type}. Valid options are 'linear' and 'non_linear'."
        #         raise Exception(err_msg)
        #     return solution_strategy

    return ROMSolver(model, custom_settings)