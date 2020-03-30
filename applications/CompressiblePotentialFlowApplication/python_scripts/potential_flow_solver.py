from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

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
            elif element_type == "perturbation_incompressible":
                self._SetUpIncompressiblePerturbationElement(formulation_settings)
            elif element_type == "perturbation_compressible":
                self._SetUpCompressiblePerturbationElement(formulation_settings)
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

    def _SetUpEmbeddedIncompressibleElement(self, formulation_settings):
        default_settings = KratosMultiphysics.Parameters(r"""{
            "element_type": "embedded_incompressible",
            "penalty_coefficient": 0.0

        }""")
        formulation_settings.ValidateAndAssignDefaults(default_settings)

        self.element_name = "EmbeddedIncompressiblePotentialFlowElement"
        self.condition_name = "PotentialWallCondition"
        self.process_info_data[KratosMultiphysics.FluidDynamicsApplication.PENALTY_COEFFICIENT] = formulation_settings["penalty_coefficient"].GetDouble()

    def _SetUpEmbeddedCompressibleElement(self, formulation_settings):
        default_settings = KratosMultiphysics.Parameters(r"""{
            "element_type": "embedded_compressible"
        }""")
        formulation_settings.ValidateAndAssignDefaults(default_settings)

        self.element_name = "EmbeddedCompressiblePotentialFlowElement"
        self.condition_name = "PotentialWallCondition"

def CreateSolver(model, custom_settings):
    return PotentialFlowSolver(model, custom_settings)

class PotentialFlowSolver(FluidSolver):

    @classmethod
    def GetDefaultSettings(cls):
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
            "maximum_iterations": 10,
            "echo_level": 0,
            "relative_tolerance": 1e-12,
            "absolute_tolerance": 1e-12,
            "compute_reactions": false,
            "reform_dofs_at_each_step": false,
            "calculate_solution_norm": false,
            "linear_solver_settings": {
                "solver_type": "amgcl"
            },
            "volume_model_part_name": "volume_model_part",
            "skin_parts":[""],
            "assign_neighbour_elements_to_conditions": false,
            "no_skin_parts": [""],
            "move_mesh_flag": false,
            "reference_chord": 1.0,
            "auxiliary_variables_list" : []
        }''')

        default_settings.AddMissingParameters(super(PotentialFlowSolver, cls).GetDefaultSettings())
        return default_settings

    def __init__(self, model, custom_settings):

        self._validate_settings_in_baseclass=True # To be removed eventually
        super(PotentialFlowSolver, self).__init__(model, custom_settings)

        # There is only a single rank in OpenMP, we always print
        self._is_printing_rank = True

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

        #construct the linear solvers
        import KratosMultiphysics.python_linear_solver_factory as linear_solver_factory
        self.linear_solver = linear_solver_factory.ConstructSolver(self.settings["linear_solver_settings"])

    def AddVariables(self):
        # Degrees of freedom
        self.main_model_part.AddNodalSolutionStepVariable(KCPFApp.VELOCITY_POTENTIAL)
        self.main_model_part.AddNodalSolutionStepVariable(KCPFApp.AUXILIARY_VELOCITY_POTENTIAL)

        # Add variables that the user defined in the ProjectParameters
        for i in range(self.settings["auxiliary_variables_list"].size()):
            variable_name = self.settings["auxiliary_variables_list"][i].GetString()
            variable = KratosMultiphysics.KratosGlobals.GetVariable(variable_name)
            self.main_model_part.AddNodalSolutionStepVariable(variable)

        KratosMultiphysics.Logger.PrintInfo("::[PotentialFlowSolver]:: ", "Variables ADDED")

    def AddDofs(self):
        KratosMultiphysics.VariableUtils().AddDof(KCPFApp.VELOCITY_POTENTIAL, self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(KCPFApp.AUXILIARY_VELOCITY_POTENTIAL, self.main_model_part)

    def Initialize(self):
        self._ComputeNodalNeighbours()

        time_scheme = KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme()
        strategy = self._GetStrategyType()
        if strategy == "linear":
            # TODO: Rename to self.strategy once we upgrade the base FluidDynamicsApplication solvers
            self.solver = KratosMultiphysics.ResidualBasedLinearStrategy(
                self.GetComputingModelPart(),
                time_scheme,
                self.linear_solver,
                self.settings["compute_reactions"].GetBool(),
                self.settings["reform_dofs_at_each_step"].GetBool(),
                self.settings["calculate_solution_norm"].GetBool(),
                self.settings["move_mesh_flag"].GetBool())
        elif strategy == "non_linear":
            conv_criteria = KratosMultiphysics.ResidualCriteria(
                self.settings["relative_tolerance"].GetDouble(),
                self.settings["absolute_tolerance"].GetDouble())
            max_iterations = self.settings["maximum_iterations"].GetInt()

            self.solver = KratosMultiphysics.ResidualBasedNewtonRaphsonStrategy(
                self.GetComputingModelPart(),
                time_scheme,
                self.linear_solver,
                conv_criteria,
                max_iterations,
                self.settings["compute_reactions"].GetBool(),
                self.settings["reform_dofs_at_each_step"].GetBool(),
                self.settings["move_mesh_flag"].GetBool())
        else:
            raise Exception("Element not implemented")

        (self.solver).SetEchoLevel(self.settings["echo_level"].GetInt())
        self.solver.Initialize()

    def AdvanceInTime(self, current_time):
        raise Exception("AdvanceInTime is not implemented. Potential Flow simulations are steady state.")

    def _ComputeNodalNeighbours(self):
        # Find nodal neigbours util call
        KratosMultiphysics.FindNodalNeighboursProcess(self.main_model_part).Execute()

    def _GetStrategyType(self):
        element_type = self.settings["formulation"]["element_type"].GetString()
        if "incompressible" in element_type:
            if not self.settings["formulation"].Has("penalty_coefficient"):
                strategy = "linear"
            elif self.settings["formulation"]["penalty_coefficient"].GetDouble() == 0.0:
                strategy = "linear"
            else:
                strategy = "non_linear"
        elif "compressible" in element_type:
            strategy = "non_linear"
        else:
            strategy = ""

        return strategy
