import KratosMultiphysics as kratos
import KratosMultiphysics.ConvectionDiffusionApplication as convdiff

from python_solver import PythonSolver

def CreateSolver(model, settings):
    return AdjointDiffusionSolver(model, settings)

class AdjointDiffusionSolver(PythonSolver):

    def __init__(self,model,settings):
        super(AdjointDiffusionSolver,self).__init__(model, settings)

        self.ValidateSettings()

        domain_size = self.settings["domain_size"].GetInt()
        if domain_size not in [2,3]:
            raise Exception("Unsupported domain_size: ", domain_size)

        model_part_name = self.settings["model_part_name"].GetString()
        if model_part_name == "":
            raise Exception("Empty model_part_name provided")

        if self.model.HasModelPart(model_part_name):
            self.model_part = self.model.GetModelPart(model_part_name)
        else:
            self.model_part = self.model.CreateModelPart(model_part_name)

        self.model_part.ProcessInfo[kratos.DOMAIN_SIZE] = domain_size

    def AddVariables(self):
        self.model_part.AddNodalSolutionStepVariable(kratos.TEMPERATURE)
        self.model_part.AddNodalSolutionStepVariable(kratos.HEAT_FLUX) # needed by primal element RHS
        self.model_part.AddNodalSolutionStepVariable(convdiff.ADJOINT_HEAT_TRANSFER)

    def AddDofs(self):
        for node in self.model_part.Nodes:
            node.AddDof(convdiff.ADJOINT_HEAT_TRANSFER)

    def GetComputingModelPart(self):
        return self.model_part

    def ImportModelPart(self):
        # we can use the default implementation in the base class
        self._ImportModelPart(self.model_part,self.settings["model_import_settings"])

    def PrepareModelPart(self):
        if not self.model_part.ProcessInfo[kratos.IS_RESTARTED]:

            throw_errors = False
            kratos.TetrahedralMeshOrientationCheck(self.model_part, throw_errors).Execute()

            kratos.ReplaceElementsAndConditionsProcess(self.model_part,self.settings["element_replace_settings"]).Execute()

    def Initialize(self):

        domain_size = self.model_part.ProcessInfo[kratos.DOMAIN_SIZE]

        if self.settings["response_function_settings"]["response_type"].GetString() == "point_temperature":
            self.response_function = convdiff.PointTemperatureResponseFunction(self.settings["response_function_settings"]["custom_settings"],self.model_part)
        else:
            raise Exception("invalid response_type: " + self.settings["response_function_settings"]["response_type"].GetString())

        self.sensitivity_builder = kratos.SensitivityBuilder(self.settings["sensitivity_settings"], self.model_part, self.response_function)


        import KratosMultiphysics.python_linear_solver_factory as linear_solver_factory
        self.linear_solver = linear_solver_factory.ConstructSolver(self.settings["linear_solver_settings"])

        self.time_scheme = kratos.ResidualBasedAdjointStaticScheme(self.response_function)

        builder_and_solver = kratos.ResidualBasedBlockBuilderAndSolver(self.linear_solver)

        self.solver = kratos.ResidualBasedLinearStrategy(self.model_part,
                                                         self.time_scheme,
                                                         self.linear_solver,
                                                         builder_and_solver,
                                                         False,
                                                         False,
                                                         False,
                                                         False)

        self.solver.SetEchoLevel(self.settings["echo_level"].GetInt())

        self.solver.Initialize()
        self.response_function.Initialize()
        self.sensitivity_builder.Initialize()
        kratos.Logger.PrintInfo(self.__class__.__name__, "Solver initialization finished.")

    def InitializeSolutionStep(self):
        self.solver.InitializeSolutionStep()
        self.response_function.InitializeSolutionStep()

    def Predict(self):
        self.solver.Predict()

    def SolveSolutionStep(self):
        return self.solver.SolveSolutionStep()

    def FinalizeSolutionStep(self):
        (self.solver).FinalizeSolutionStep()
        self.response_function.FinalizeSolutionStep()
        self.sensitivity_builder.UpdateSensitivities()

    def Check(self):
        (self.solver).Check()

    def Clear(self):
        (self.solver).Clear()

    def GetDefaultSettings(self):
        return kratos.Parameters(r'''{
            "solver_type" : "adjoint_stationary",
            "model_part_name": "",
            "domain_size": 0,
            "response_function_settings" : {
                "response_type" : "point_temperature"
            },
            "sensitivity_settings" : {},
            "model_import_settings" : {
                "input_type"     : "mdpa",
                "input_filename" : ""
            },
            "linear_solver_settings" : {
                "solver_type" : "amgcl"
            },
            "element_replace_settings" : {
                "element_name" : "AdjointHeatDiffusionElement2D3N",
                "condition_name" : ""
            },
            "volume_model_part_name" : "volume_model_part",
            "domain_model_parts"                 : [],
            "boundary_model_parts"               : [],
            "echo_level"  : 0
        }''')
