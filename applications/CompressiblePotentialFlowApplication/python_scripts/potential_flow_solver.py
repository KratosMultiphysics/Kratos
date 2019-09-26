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

        if formulation_settings.Has("element_type"):
            element_type = formulation_settings["element_type"].GetString()
            if element_type == "incompressible":
                self._SetUpIncompressibleElement(formulation_settings)
            elif element_type == "compressible":
                self._SetUpCompressibleElement(formulation_settings)
            elif element_type == "embedded_incompressible":
                self._SetUpEmbeddedIncompressibleElement(formulation_settings)
        else:
            raise RuntimeError("Argument \'element_type\' not found in formulation settings.")

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

    def _SetUpEmbeddedIncompressibleElement(self, formulation_settings):
        default_settings = KratosMultiphysics.Parameters(r"""{
            "element_type": "embedded_incompressible"
        }""")
        formulation_settings.ValidateAndAssignDefaults(default_settings)

        self.element_name = "EmbeddedIncompressiblePotentialFlowElement"
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
                "materials_filename": "unknown_materials.json"
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
        self.min_buffer_size = 1
        self.domain_size = custom_settings["domain_size"].GetInt()
        self.reference_chord = custom_settings["reference_chord"].GetDouble()
        self.main_model_part.ProcessInfo.SetValue(KCPFApp.REFERENCE_CHORD,self.reference_chord)

        #construct the linear solvers
        import KratosMultiphysics.python_linear_solver_factory as linear_solver_factory
        self.linear_solver = linear_solver_factory.ConstructSolver(self.settings["linear_solver_settings"])

    def AddVariables(self):
        # Degrees of freedom
        self.main_model_part.AddNodalSolutionStepVariable(KCPFApp.VELOCITY_POTENTIAL)
        self.main_model_part.AddNodalSolutionStepVariable(KCPFApp.AUXILIARY_VELOCITY_POTENTIAL)
        self.main_model_part.AddNodalSolutionStepVariable(KCPFApp.PSI)

        # Add variables that the user defined in the ProjectParameters
        for i in range(self.settings["auxiliary_variables_list"].size()):
            variable_name = self.settings["auxiliary_variables_list"][i].GetString()
            variable = KratosMultiphysics.KratosGlobals.GetVariable(variable_name)
            self.main_model_part.AddNodalSolutionStepVariable(variable)

        KratosMultiphysics.Logger.PrintInfo("::[PotentialFlowSolver]:: ", "Variables ADDED")

    def AddDofs(self):
        KratosMultiphysics.VariableUtils().AddDof(KCPFApp.VELOCITY_POTENTIAL, self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(KCPFApp.AUXILIARY_VELOCITY_POTENTIAL, self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(KCPFApp.PSI, self.main_model_part)

    def Initialize(self):
        self._ComputeNodalNeighbours()

        time_scheme = KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme()
        if "incompressible" in self.settings["formulation"]["element_type"].GetString():
            # TODO: Rename to self.strategy once we upgrade the base FluidDynamicsApplication solvers
            self.solver = KratosMultiphysics.ResidualBasedLinearStrategy(
                self.GetComputingModelPart(),
                time_scheme,
                self.linear_solver,
                self.settings["compute_reactions"].GetBool(),
                self.settings["reform_dofs_at_each_step"].GetBool(),
                self.settings["calculate_solution_norm"].GetBool(),
                self.settings["move_mesh_flag"].GetBool())
        elif "compressible" in self.settings["formulation"]["element_type"].GetString():
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
        avg_elem_num = 10
        avg_node_num = 10
        KratosMultiphysics.FindNodalNeighboursProcess(
            self.main_model_part, avg_elem_num, avg_node_num).Execute()

    def SolveSolutionStep(self):
        import time as time
        start_time = time.time()
        super(PotentialFlowSolver, self).SolveSolutionStep()
        exe_time = time.time() - start_time
        print('Executing SolveSolutionStep took ' + str(round(exe_time/60, 2)) + ' min')
        # start_time = time.time()
        # self._ComputeConditionNumber()
        # exe_time = time.time() - start_time
        # print('Executing _ComputeConditionNumber took ' + str(round(exe_time, 2)) + ' sec')

    def _ComputeConditionNumber(self):
        KratosMultiphysics.Logger.PrintInfo("::[PotentialFlowSolver]:: ", "Computing Condition Number")

        import eigen_solver_factory
        settings_max = KratosMultiphysics.Parameters("""
        {
            "solver_type"             : "power_iteration_highest_eigenvalue_solver",
            "max_iteration"           : 10000,
            "tolerance"               : 1e-9,
            "required_eigen_number"   : 1,
            "verbosity"               : 0,
            "linear_solver_settings"  : {
                "solver_type"             : "ExternalSolversApplication.super_lu",
                "max_iteration"           : 500,
                "tolerance"               : 1e-9,
                "scaling"                 : false,
                "verbosity"               : 0
            }
        }
        """)
        eigen_solver_max = eigen_solver_factory.ConstructSolver(settings_max)
        settings_min = KratosMultiphysics.Parameters("""
        {
            "solver_type"             : "power_iteration_eigenvalue_solver",
            "max_iteration"           : 10000,
            "tolerance"               : 1e-9,
            "required_eigen_number"   : 1,
            "verbosity"               : 0,
            "linear_solver_settings"  : {
                "solver_type"             : "ExternalSolversApplication.super_lu",
                "max_iteration"           : 500,
                "tolerance"               : 1e-9,
                "scaling"                 : false,
                "verbosity"               : 0
            }
        }
        """)

        eigen_solver_min = eigen_solver_factory.ConstructSolver(settings_min)
        condition_number = KratosMultiphysics.ConditionNumberUtility().GetConditionNumber(self.solver.GetSystemMatrix(), eigen_solver_max, eigen_solver_min)
        KratosMultiphysics.Logger.PrintInfo(' condition_number = ', "{:.2e}".format(condition_number))