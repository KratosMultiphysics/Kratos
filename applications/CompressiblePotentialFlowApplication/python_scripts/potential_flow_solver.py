# makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division

# Importing the Kratos Library
import KratosMultiphysics
from python_solver import PythonSolver

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
            "model_part_name"        : "model",
            "domain_size"            : 2,
            "solver_type": "potential_flow_solver",
            "domain_size"     : 2,
            "model_part_name" : "MainModelPart",
            "echo_level": 1,
            "relative_tolerance": 1e-5,
            "absolute_tolerance": 1e-9,
            "maximum_iterations": 1,
            "compute_reactions": false,
            "reform_dofs_at_each_step": false,
            "calculate_solution_norm" : false,
            "volume_model_part_name" : "volume_model_part",
            "skin_parts":[],
            "no_skin_parts"                : [],
            "dimension"             : 0,
            "node_id"               : 0,
            "epsilon"               : 1e-6,
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
                    "gmres_krylov_space_dimension": 500,
                    "smoother_type":"ilu0",
                    "coarsening_type":"ruge_stuben",
                    "coarse_enough" : 5000,
                    "krylov_type": "lgmres",
                    "tolerance": 1e-9,
                    "verbosity": 0,
                    "scaling": false
            }


        }""")

        # Overwrite the default settings with user-provided parameters.
        self.settings = custom_settings
        self.settings.ValidateAndAssignDefaults(default_settings)

        model_part_name = self.settings["model_part_name"].GetString()

        if model_part_name == "":
            raise Exception('Please provide the model part name as the "model_part_name" (string) parameter!')

        if self.model.HasModelPart(model_part_name):
            self.main_model_part = self.model.GetModelPart(model_part_name)
            self.solver_imports_model_part = False
        else:
            self.main_model_part = model.CreateModelPart(model_part_name)
            self.solver_imports_model_part = True

        self.domain_size = custom_settings["domain_size"].GetInt()
        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, self.domain_size)
        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DENSITY, 1.225)
        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.MIU,5)#geometry angle
        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.LAMBDA, 1.4)
        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.SOUND_VELOCITY, 340.0)


        #construct the linear solvers
        import KratosMultiphysics.python_linear_solver_factory as linear_solver_factory
        self.linear_solver = linear_solver_factory.ConstructSolver(self.settings["linear_solver_settings"])
        self.perturbate_model_part=False

    def AddVariables(self):
        # Degrees of freedom
        self.main_model_part.AddNodalSolutionStepVariable(KCPFApp.VELOCITY_POTENTIAL)
        self.main_model_part.AddNodalSolutionStepVariable(KCPFApp.AUXILIARY_VELOCITY_POTENTIAL)

        # Kratos variables
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.TEMPERATURE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL)

    def AddDofs(self):
        KratosMultiphysics.VariableUtils().AddDof(KCPFApp.VELOCITY_POTENTIAL, self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(KCPFApp.AUXILIARY_VELOCITY_POTENTIAL, self.main_model_part)

    def Initialize(self):

        time_scheme = KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme()


        if self.settings["element_replace_settings"]["element_name"].GetString() == "CompressiblePotentialFlowElement2D3N":
            conv_criteria = KratosMultiphysics.ResidualCriteria(
                self.settings["relative_tolerance"].GetDouble(),
                self.settings["absolute_tolerance"].GetDouble())
            max_iterations = self.settings["maximum_iterations"].GetInt()

            self.solver = KratosMultiphysics.ResidualBasedNewtonRaphsonStrategy(
                self.main_model_part,
                time_scheme,
                self.linear_solver,
                conv_criteria,
                max_iterations,
                self.settings["compute_reactions"].GetBool(),
                self.settings["reform_dofs_at_each_step"].GetBool(),
                self.move_mesh_flag)
        else:
            builder_and_solver = KratosMultiphysics.ResidualBasedBlockBuilderAndSolver(self.linear_solver)
            self.solver = KratosMultiphysics.ResidualBasedLinearStrategy(
                self.main_model_part,
                time_scheme,
                self.linear_solver,
                builder_and_solver,
                self.settings["compute_reactions"].GetBool(),
                self.settings["reform_dofs_at_each_step"].GetBool(),
                self.settings["calculate_solution_norm"].GetBool(),
                self.move_mesh_flag)

        (self.solver).SetEchoLevel(self.settings["echo_level"].GetInt())
        self.solver.Check()
        self.solver.Initialize()
    def PrepareModelPart(self):
        if not self.model.HasModelPart(self.settings["model_part_name"].GetString()):
            self.model.AddModelPart(self.main_model_part)


    def ImportModelPart(self):

        """This function imports the ModelPart
        """
        if self.solver_imports_model_part:
            if(self.settings["model_import_settings"]["input_type"].GetString() == "mdpa"):
                #here it would be the place to import restart data if required

                self._ImportModelPart(self.main_model_part,self.settings["model_import_settings"])

                # IOdir=self.settings["model_import_settings"]["input_filename"].GetString()
                # KratosMultiphysics.ModelPartIO(IOdir).ReadModelPart(self.main_model_part)

                throw_errors = False
                KratosMultiphysics.TetrahedralMeshOrientationCheck(self.main_model_part,throw_errors).Execute()
                #here we replace the dummy elements we read with proper elements
            
            KratosMultiphysics.ReplaceElementsAndConditionsProcess(self.main_model_part, self.settings["element_replace_settings"]).Execute()

        # else:
        # raise Exception("other input options are not yet implemented")
        current_buffer_size = self.main_model_part.GetBufferSize()
        if(self.GetMinimumBufferSize() > current_buffer_size):
            self.main_model_part.SetBufferSize( self.GetMinimumBufferSize() )

        print ("model reading finished")

    def GetMinimumBufferSize(self):
        return 1

    def GetComputingModelPart(self):
        return self.main_model_part

    def InitializeSolutionStep(self):
        self.solver.InitializeSolutionStep()

    def Predict(self):
        self.solver.Predict()

    def SolveSolutionStep(self):
        self.solver.SolveSolutionStep()

    def FinalizeSolutionStep(self):
        self.solver.FinalizeSolutionStep()

    def SetEchoLevel(self, level):
        self.solver.SetEchoLevel(level)

    def Clear(self):
        self.solver.Clear()

    def AdvanceInTime(self, current_time):
        raise Exception("AdvanceInTime is not implemented. Potential Flow simulations are steady state.")

