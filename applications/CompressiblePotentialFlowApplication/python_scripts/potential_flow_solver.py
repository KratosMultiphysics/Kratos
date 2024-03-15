# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.CompressiblePotentialFlowApplication as KCPFApp

# Importing the base class
from KratosMultiphysics.FluidDynamicsApplication.fluid_solver import FluidSolver

class PotentialFlowFormulation(object):
    def __init__(self, formulation_settings):
        self.element_name = None
        self.condition_name = None
        self.process_info_data = {}

        if formulation_settings.Has("element_type"):
            element_type = formulation_settings["element_type"].GetString()
            if element_type == "incompressible":
                self._SetUpIncompressibleElement(formulation_settings)
            elif element_type == "compressible":
                self._SetUpCompressibleElement(formulation_settings)
            elif element_type == "embedded_incompressible":
                self._SetUpEmbeddedIncompressibleElement(formulation_settings)
            elif element_type == "embedded_compressible":
                self._SetUpEmbeddedCompressibleElement(formulation_settings)
            elif element_type == "embedded_perturbation_transonic":
                self._SetUpEmbeddedTransonicPerturbationElement(formulation_settings)
            elif element_type == "perturbation_incompressible":
                self._SetUpIncompressiblePerturbationElement(formulation_settings)
            elif element_type == "perturbation_compressible":
                self._SetUpCompressiblePerturbationElement(formulation_settings)
            elif element_type == "perturbation_transonic":
                self._SetUpTransonicPerturbationElement(formulation_settings)
        else:
            raise RuntimeError("Argument \'element_type\' not found in formulation settings.")

    def SetProcessInfo(self, model_part):
        for variable,value in self.process_info_data.items():
            model_part.ProcessInfo[variable] = value

    def _SetUpIncompressibleElement(self, formulation_settings):
        default_settings = KratosMultiphysics.Parameters(r"""{
            "element_type": "incompressible"
        }""")
        formulation_settings.ValidateAndAssignDefaults(default_settings)

        self.element_name = "IncompressiblePotentialFlowElement"
        self.condition_name = "PotentialWallCondition"

    def _SetUpCompressibleElement(self, formulation_settings):
        default_settings = KratosMultiphysics.Parameters(r"""{
            "element_type": "compressible"
        }""")
        formulation_settings.ValidateAndAssignDefaults(default_settings)

        self.element_name = "CompressiblePotentialFlowElement"
        self.condition_name = "PotentialWallCondition"

    def _SetUpIncompressiblePerturbationElement(self, formulation_settings):
        default_settings = KratosMultiphysics.Parameters(r"""{
            "element_type": "perturbation_incompressible"
        }""")
        formulation_settings.ValidateAndAssignDefaults(default_settings)

        self.element_name = "IncompressiblePerturbationPotentialFlowElement"
        self.condition_name = "PotentialWallCondition"

    def _SetUpCompressiblePerturbationElement(self, formulation_settings):
        default_settings = KratosMultiphysics.Parameters(r"""{
            "element_type": "perturbation_compressible"
        }""")
        formulation_settings.ValidateAndAssignDefaults(default_settings)

        self.element_name = "CompressiblePerturbationPotentialFlowElement"
        self.condition_name = "PotentialWallCondition"

    def _SetUpTransonicPerturbationElement(self, formulation_settings):
        default_settings = KratosMultiphysics.Parameters(r"""{
            "element_type": "perturbation_transonic"
        }""")
        formulation_settings.ValidateAndAssignDefaults(default_settings)

        self.element_name = "TransonicPerturbationPotentialFlowElement"
        self.condition_name = "PotentialWallCondition"

    def _SetUpEmbeddedIncompressibleElement(self, formulation_settings):
        default_settings = KratosMultiphysics.Parameters(r"""{
            "element_type": "embedded_incompressible",
            "stabilization_factor": 0.0,
            "penalty_coefficient": 0.0

        }""")
        formulation_settings.ValidateAndAssignDefaults(default_settings)

        self.element_name = "EmbeddedIncompressiblePotentialFlowElement"
        self.condition_name = "PotentialWallCondition"
        self.process_info_data[KratosMultiphysics.STABILIZATION_FACTOR] = formulation_settings["stabilization_factor"].GetDouble()
        self.process_info_data[KratosMultiphysics.PENALTY_COEFFICIENT] = formulation_settings["penalty_coefficient"].GetDouble()

    def _SetUpEmbeddedCompressibleElement(self, formulation_settings):
        default_settings = KratosMultiphysics.Parameters(r"""{
            "element_type": "embedded_compressible",
            "stabilization_factor": 0.0,
            "penalty_coefficient": 0.0
        }""")
        formulation_settings.ValidateAndAssignDefaults(default_settings)

        self.element_name = "EmbeddedCompressiblePotentialFlowElement"
        self.condition_name = "PotentialWallCondition"
        self.process_info_data[KratosMultiphysics.STABILIZATION_FACTOR] = formulation_settings["stabilization_factor"].GetDouble()
        self.process_info_data[KratosMultiphysics.PENALTY_COEFFICIENT] = formulation_settings["penalty_coefficient"].GetDouble()

    def _SetUpEmbeddedTransonicPerturbationElement(self, formulation_settings):
        default_settings = KratosMultiphysics.Parameters(r"""{
            "element_type": "embedded_perturbation_transonic",
            "stabilization_factor": 0.0,
            "penalty_coefficient": 0.0
        }""")
        formulation_settings.ValidateAndAssignDefaults(default_settings)

        self.element_name = "EmbeddedTransonicPerturbationPotentialFlowElement"
        self.condition_name = "PotentialWallCondition"
        self.process_info_data[KratosMultiphysics.STABILIZATION_FACTOR] = formulation_settings["stabilization_factor"].GetDouble()
        self.process_info_data[KratosMultiphysics.PENALTY_COEFFICIENT] = formulation_settings["penalty_coefficient"].GetDouble()


def CreateSolver(model, custom_settings):
    return PotentialFlowSolver(model, custom_settings)

class PotentialFlowSolver(FluidSolver):

    @classmethod
    def GetDefaultParameters(cls):
        # Default settings string in json format
        default_settings = KratosMultiphysics.Parameters(r'''{
            "solver_type": "potential_flow_solver",
            "model_part_name": "PotentialFluidModelPart",
            "domain_size": -1,
            "model_import_settings": {
                "input_type": "mdpa",
                "input_filename": "unknown_name"
            },
            "material_import_settings": {
                "materials_filename": ""
            },
            "formulation": {
                "element_type": "incompressible"
            },
            "element_replace_settings": {
                "condition_name":  "",
                "element_name": ""
            },
            "maximum_iterations": 10,
            "echo_level": 0,
            "potential_application_echo_level": 0,
            "convergence_criterion": "residual_criterion",
            "solving_strategy_settings": {
                "type" : "newton_raphson",
                "advanced_settings" : {}
            },
            "relative_tolerance": 1e-12,
            "absolute_tolerance": 1e-12,
            "compute_reactions": false,
            "reform_dofs_at_each_step": false,
            "linear_solver_settings": {
                "solver_type": "amgcl"
            },
            "volume_model_part_name": "volume_model_part",
            "skin_parts":[""],
            "assign_neighbour_elements_to_conditions": false,
            "no_skin_parts": [""],
            "time_stepping"                : {
                "automatic_time_step" : false,
                "CFL_number"          : 1,
                "minimum_delta_time"  : 1e-4,
                "maximum_delta_time"  : 1.0,
                "time_step":            1.0
            },
            "move_mesh_flag": false,
            "reference_chord": 1.0,
            "auxiliary_variables_list" : []
        }''')

        default_settings.AddMissingParameters(super(PotentialFlowSolver, cls).GetDefaultParameters())
        return default_settings

    def __init__(self, model, custom_settings):

        self._validate_settings_in_baseclass=True # To be removed eventually
        super().__init__(model, custom_settings)
        self._enforce_element_and_condition_replacement = True #TODO: Remove once we remove the I/O from the solver

        # Set the element and condition names for the replace settings
        self.formulation = PotentialFlowFormulation(self.settings["formulation"])
        self.element_name = self.formulation.element_name
        self.condition_name = self.formulation.condition_name
        self.formulation.SetProcessInfo(self.main_model_part)
        self.min_buffer_size = 1
        self.domain_size = custom_settings["domain_size"].GetInt()
        self.reference_chord = custom_settings["reference_chord"].GetDouble()
        self.main_model_part.ProcessInfo.SetValue(KCPFApp.REFERENCE_CHORD,self.reference_chord)
        self.element_has_nodal_properties = False

    def AddVariables(self):
        # Degrees of freedom
        self.main_model_part.AddNodalSolutionStepVariable(KCPFApp.VELOCITY_POTENTIAL)
        self.main_model_part.AddNodalSolutionStepVariable(KCPFApp.AUXILIARY_VELOCITY_POTENTIAL)
        self.main_model_part.AddNodalSolutionStepVariable(KCPFApp.REACTION_VELOCITY_POTENTIAL)
        self.main_model_part.AddNodalSolutionStepVariable(KCPFApp.REACTION_AUXILIARY_VELOCITY_POTENTIAL)

        # Add variables that the user defined in the ProjectParameters
        for i in range(self.settings["auxiliary_variables_list"].size()):
            variable_name = self.settings["auxiliary_variables_list"][i].GetString()
            variable = KratosMultiphysics.KratosGlobals.GetVariable(variable_name)
            self.main_model_part.AddNodalSolutionStepVariable(variable)

        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Variables ADDED")

    def AddDofs(self):
        dofs_and_reactions_to_add = []
        dofs_and_reactions_to_add.append(["VELOCITY_POTENTIAL", "REACTION_VELOCITY_POTENTIAL"])
        dofs_and_reactions_to_add.append(["AUXILIARY_VELOCITY_POTENTIAL", "REACTION_AUXILIARY_VELOCITY_POTENTIAL"])
        KratosMultiphysics.VariableUtils.AddDofsList(dofs_and_reactions_to_add, self.main_model_part)

    def Initialize(self):
        self._ComputeNodalElementalNeighbours()

        solution_strategy = self._GetSolutionStrategy()
        solution_strategy.SetEchoLevel(self.settings["echo_level"].GetInt())
        solution_strategy.Initialize()
        self.GetComputingModelPart().ProcessInfo.SetValue(
            KCPFApp.ECHO_LEVEL, self.settings["potential_application_echo_level"].GetInt())

        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Solver initialization finished.")

    def _ComputeNodalElementalNeighbours(self):
        # Find nodal neigbours util call
        KratosMultiphysics.FindGlobalNodalElementalNeighboursProcess(self.GetComputingModelPart()).Execute()

    def _GetStrategyType(self):
        element_type = self.settings["formulation"]["element_type"].GetString()
        if "incompressible" in element_type:
            if not self.settings["formulation"].Has("stabilization_factor"):
                strategy_type = "linear"
            elif self.settings["formulation"]["stabilization_factor"].GetDouble() == 0.0:
                strategy_type = "linear"
            else:
                strategy_type = "non_linear"
        elif "compressible" or "transonic" in element_type:
            strategy_type = "non_linear"
        else:
            strategy_type = None
        return strategy_type

    def _CreateBuilderAndSolver(self):
        linear_solver = self._GetLinearSolver()
        return KratosMultiphysics.ResidualBasedBlockBuilderAndSolver(linear_solver)

    def _CreateScheme(self):
        # Fake scheme creation to do the solution update
        scheme = KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme()
        return scheme

    def _CreateConvergenceCriterion(self):
        criterion = self.settings["convergence_criterion"].GetString()
        if criterion == "solution_criterion":
            convergence_criterion = KratosMultiphysics.DisplacementCriteria(
                self.settings["relative_tolerance"].GetDouble(),
                self.settings["absolute_tolerance"].GetDouble())
        elif criterion == "residual_criterion":
            convergence_criterion = KratosMultiphysics.ResidualCriteria(
                self.settings["relative_tolerance"].GetDouble(),
                self.settings["absolute_tolerance"].GetDouble())
        else:
            err_msg =  "The requested convergence criterion \"" + criterion + "\" is not available!\n"
            err_msg += "Available options are: \"solution_criterion\", \"residual_criterion\""
            raise Exception(err_msg)
        convergence_criterion.SetEchoLevel(self.settings["echo_level"].GetInt())
        return convergence_criterion

    def _CreateSolutionStrategy(self):
        strategy_type = self._GetStrategyType()
        if strategy_type == "linear":
            solution_strategy = self._CreateLinearStrategy()
        elif strategy_type == "non_linear":
            # Create strategy
            if self.settings["solving_strategy_settings"]["type"].GetString() == "newton_raphson":
                solution_strategy = self._CreateNewtonRaphsonStrategy()
            elif self.settings["solving_strategy_settings"]["type"].GetString() == "line_search":
                solution_strategy = self._CreateLineSearchStrategy()
        else:
            err_msg = "Unknown strategy type: \'" + strategy_type + "\'. Valid options are \'linear\' and \'non_linear\'."
            raise Exception(err_msg)
        return solution_strategy

    def _CreateLineSearchStrategy(self):
        if self.settings["solving_strategy_settings"].Has("advanced_settings"):
            settings = self.settings["solving_strategy_settings"]["advanced_settings"]
            settings.AddMissingParameters(self._GetDefaultLineSearchParameters())
        else:
            settings = self._GetDefaultLineSearchParameters()
        settings.AddValue("max_iteration", self.settings["maximum_iterations"])
        settings.AddValue("compute_reactions", self.settings["compute_reactions"])
        settings.AddValue("reform_dofs_at_each_step", self.settings["reform_dofs_at_each_step"])
        settings.AddValue("move_mesh_flag", self.settings["move_mesh_flag"])
        computing_model_part = self.GetComputingModelPart()
        time_scheme = self._GetScheme()
        convergence_criterion = self._GetConvergenceCriterion()
        builder_and_solver = self._GetBuilderAndSolver()
        solution_strategy = KratosMultiphysics.LineSearchStrategy(computing_model_part,
            time_scheme,
            convergence_criterion,
            builder_and_solver,
            settings)
        return solution_strategy

    @classmethod
    def _GetDefaultLineSearchParameters(self):
        default_line_search_parameters = KratosMultiphysics.Parameters(r"""{
                "max_line_search_iterations" : 5,
                "first_alpha_value"          : 0.5,
                "second_alpha_value"         : 1.0,
                "min_alpha"                  : 0.1,
                "max_alpha"                  : 2.0,
                "line_search_tolerance"      : 0.5
            }""")
        return default_line_search_parameters

    def _SetPhysicalProperties(self):
        # There are no properties in the potential flow solver. Free stream quantities are defined in the apply_far_field_process.py
        return True
