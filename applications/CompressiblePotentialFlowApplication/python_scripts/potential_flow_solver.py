# makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division

# Importing the Kratos Library
import KratosMultiphysics

import KratosMultiphysics.CompressiblePotentialFlowApplication as KCPFApp

# Importing the base class
from KratosMultiphysics.python_solver import PythonSolver


def CreateSolver(model, custom_settings):
    return LaplacianSolver(model, custom_settings)


class LaplacianSolver(PythonSolver):
    def __init__(self, model, custom_settings):
        super(LaplacianSolver, self).__init__(model, custom_settings)
        self.move_mesh_flag = False

        ##settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
        {
            "solver_type": "potential_flow_solver",
            "domain_size"     : 2,
            "model_part_name" : "MainModelPart",
            "echo_level": 1,
            "relative_tolerance": 1e-5,
            "absolute_tolerance": 1e-9,
            "maximum_iterations": 10,
            "compute_reactions": false,
            "reform_dofs_at_each_step": false,
            "calculate_solution_norm" : false,
            "volume_model_part_name" : "volume_model_part",
            "skin_parts":[],
            "no_skin_parts"                : [],
            "model_import_settings": {
                    "input_type": "mdpa",
                    "input_filename": "unknown_name"
            },
            "element_replace_settings": {
                    "element_name":"IncompressiblePotentialFlowElement2D3N",
                    "condition_name": "PotentialWallCondition2D2N"
            },
            "linear_solver_settings": {
                    "solver_type": "amgcl",
                    "max_iteration": 400,
                    "gmres_krylov_space_dimension": 100,
                    "smoother_type":"ilu0",
                    "coarsening_type":"ruge_stuben",
                    "coarse_enough" : 5000,
                    "krylov_type": "lgmres",
                    "tolerance": 1e-9,
                    "verbosity": 3,
                    "scaling": false
            }
        }""")

        # Overwrite the default settings with user-provided parameters.
        self.settings = custom_settings
        self.settings.ValidateAndAssignDefaults(default_settings)
        model_part_name = self.settings["model_part_name"].GetString()

        if model_part_name == "":
            raise Exception('Please specify a model_part name!')

        # This will be changed once the Model is fully supported!
        if self.model.HasModelPart(model_part_name):
            self.main_model_part = self.model[model_part_name]
        else:
            self.main_model_part = self.model.CreateModelPart(model_part_name)
            domain_size = self.settings["domain_size"].GetInt()
            if domain_size < 0:
                raise Exception('Please specify a "domain_size" >= 0!')
            self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, domain_size)

        # Construct the linear solvers
        import KratosMultiphysics.python_linear_solver_factory as linear_solver_factory
        self.linear_solver = linear_solver_factory.ConstructSolver(self.settings["linear_solver_settings"])

        print("Construction of LaplacianSolver finished")

    def AddVariables(self):
        # Degrees of freedom
        self.main_model_part.AddNodalSolutionStepVariable(KCPFApp.VELOCITY_POTENTIAL)
        self.main_model_part.AddNodalSolutionStepVariable(KCPFApp.AUXILIARY_VELOCITY_POTENTIAL)

        # Kratos variables
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL)

    def AddDofs(self):
        KratosMultiphysics.VariableUtils().AddDof(KCPFApp.VELOCITY_POTENTIAL, self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(KCPFApp.AUXILIARY_VELOCITY_POTENTIAL, self.main_model_part)

    def Initialize(self):

        time_scheme = KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme()

        if(self.settings["element_replace_settings"]["element_name"].GetString() == "IncompressiblePotentialFlowElement2D3N"):
            self.solution_strategy = KratosMultiphysics.ResidualBasedLinearStrategy(
                self.main_model_part,
                time_scheme,
                self.linear_solver,
                self.settings["compute_reactions"].GetBool(),
                self.settings["reform_dofs_at_each_step"].GetBool(),
                self.settings["calculate_solution_norm"].GetBool(),
                self.move_mesh_flag)
        elif(self.settings["element_replace_settings"]["element_name"].GetString() == "CompressiblePotentialFlowElement2D3N"):
            conv_criteria = KratosMultiphysics.ResidualCriteria(
                self.settings["relative_tolerance"].GetDouble(), 
                self.settings["absolute_tolerance"].GetDouble())
            max_iterations = self.settings["maximum_iterations"].GetInt()

            self.solution_strategy = KratosMultiphysics.ResidualBasedNewtonRaphsonStrategy(
                self.main_model_part,
                time_scheme,
                self.linear_solver,
                conv_criteria,
                max_iterations,
                self.settings["compute_reactions"].GetBool(),
                self.settings["reform_dofs_at_each_step"].GetBool(),
                self.move_mesh_flag)
        else:
            raise Exception("Element not implemented")

        (self.solution_strategy).SetEchoLevel(self.settings["echo_level"].GetInt())
        self.solution_strategy.Initialize()

    def Check(self):
        self.solution_strategy.Check()

    def ImportModelPart(self):

        if(self.settings["model_import_settings"]["input_type"].GetString() == "mdpa"):
            #here it would be the place to import restart data if required
            print(self.settings["model_import_settings"]["input_filename"].GetString())
            KratosMultiphysics.ModelPartIO(self.settings["model_import_settings"]["input_filename"].GetString()).ReadModelPart(self.main_model_part)

            throw_errors = False
            KratosMultiphysics.TetrahedralMeshOrientationCheck(self.main_model_part,throw_errors).Execute()
            #here we replace the dummy elements we read with proper elements
            self.settings.AddEmptyValue("element_replace_settings")
            if(self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 3):
                self.settings["element_replace_settings"] = KratosMultiphysics.Parameters("""
                    {
                    "element_name":"CompressiblePotentialFlowElement3D4N",
                    "condition_name": "PotentialWallCondition3D3N"
                    }
                    """)
            elif(self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] != 2):
                raise Exception("Domain size is not 2 or 3!!")

            KratosMultiphysics.ReplaceElementsAndConditionsProcess(self.main_model_part, self.settings["element_replace_settings"]).Execute()

        else:
            raise Exception("other input options are not yet implemented")

        current_buffer_size = self.main_model_part.GetBufferSize()
        if(self.GetMinimumBufferSize() > current_buffer_size):
            self.main_model_part.SetBufferSize( self.GetMinimumBufferSize() )

        print ("model reading finished")

    def GetMinimumBufferSize(self):
        return 1

    def GetComputingModelPart(self):
        return self.main_model_part

    def InitializeSolutionStep(self):
        self.solution_strategy.InitializeSolutionStep()

    def Predict(self):
        self.solution_strategy.Predict()

    def SolveSolutionStep(self):
        self.solution_strategy.SolveSolutionStep()

    def FinalizeSolutionStep(self):
        self.solution_strategy.FinalizeSolutionStep()

    def SetEchoLevel(self, level):
        self.solution_strategy.SetEchoLevel(level)

    def Clear(self):
        self.solution_strategy.Clear()

    def AdvanceInTime(self, current_time):
        raise Exception("AdvanceInTime is not implemented. Potential Flow simulations are steady state.")

