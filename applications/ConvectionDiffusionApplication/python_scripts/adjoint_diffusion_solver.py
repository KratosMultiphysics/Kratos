import KratosMultiphysics as kratos
import KratosMultiphysics.ConvectionDiffusionApplication as convdiff

from KratosMultiphysics.python_solver import PythonSolver

def CreateSolver(model, settings):
    return AdjointDiffusionSolver(model, settings)

class AdjointDiffusionSolver(PythonSolver):

    def __init__(self, model, custom_settings):
        super().__init__(model,custom_settings)

        self.min_buffer_size = 1

        model_part_name = self.settings["model_part_name"].GetString()
        if model_part_name == "":
            raise Exception("Empty model_part_name provided")

        if self.model.HasModelPart(model_part_name):
            self.model_part = self.model.GetModelPart(model_part_name)
            self.solver_imports_model_part = False
        else:
            self.model_part = self.model.CreateModelPart(model_part_name)

            domain_size = self.settings["domain_size"].GetInt()
            if domain_size not in (2,3):
                raise Exception("Unsupported domain_size: ", domain_size)

            self.model_part.ProcessInfo[kratos.DOMAIN_SIZE] = domain_size
            self.solver_imports_model_part = True

            self.DefineConvectionDiffusionSettings(self.settings["convection_diffusion_variables"])

        self.primal_model_part_name = self.settings["primal_model_part_name"].GetString()
        if self.primal_model_part_name == "":
            raise Exception("No primal_model_part_name provided")

    @classmethod
    def GetDefaultParameters(cls):

        default_settings = kratos.Parameters(r'''{
            "solver_type" : "adjoint_stationary",
            "model_part_name": "",
            "primal_model_part_name" : "",
            "domain_size": 0,
            "model_import_settings" : {
                "input_type"     : "mdpa",
                "input_filename" : ""
            },
            "convection_diffusion_variables" : {
                "diffusion_variable"            : "CONDUCTIVITY",
                "unknown_variable"              : "TEMPERATURE",
                "volume_source_variable"        : "HEAT_FLUX",
                "surface_source_variable"       : "FACE_HEAT_FLUX"
            },
            "material_import_settings" :{
                "materials_filename": ""
            },
            "linear_solver_settings" : {
                "solver_type" : "amgcl"
            },
            "response_function_settings" : {
                "response_type" : "point_temperature"
            },
            "sensitivity_settings" : {},
            "element_replace_settings" : {
                "element_name" : "AdjointDiffusionElement",
                "condition_name" : "AdjointThermalFace"
            },
            "time_stepping" : {
                "time_step" : 0.0
            },
            "time_integration_method": "implicit"
        }''')

        default_settings.AddMissingParameters(super().GetDefaultParameters())
        return default_settings

    def AddVariables(self):
        convection_diffusion_settings = self.model_part.ProcessInfo[kratos.CONVECTION_DIFFUSION_SETTINGS]

        self.model_part.AddNodalSolutionStepVariable(convection_diffusion_settings.GetUnknownVariable())
        self.model_part.AddNodalSolutionStepVariable(convection_diffusion_settings.GetDiffusionVariable())
        self.model_part.AddNodalSolutionStepVariable(convection_diffusion_settings.GetVolumeSourceVariable())
        self.model_part.AddNodalSolutionStepVariable(convection_diffusion_settings.GetSurfaceSourceVariable())
        self.model_part.AddNodalSolutionStepVariable(convdiff.ADJOINT_HEAT_TRANSFER)
        self.model_part.AddNodalSolutionStepVariable(kratos.SHAPE_SENSITIVITY)

    def AddDofs(self):
        variable_utils = kratos.VariableUtils()
        variable_utils.AddDof(convdiff.ADJOINT_HEAT_TRANSFER, self.model_part)

    def ImportModelPart(self):
        # we can use the default implementation in the base class
        if self.solver_imports_model_part:
            self._ImportModelPart(self.model_part,self.settings["model_import_settings"])

    def PrepareModelPart(self):
        if self.solver_imports_model_part:
            # ensure that the element type is the correct one
            self._set_elements_and_conditions()

            # check mesh orientation (tetrahedral mesh orientation check)
            throw_errors = False
            kratos.TetrahedralMeshOrientationCheck(self.model_part, throw_errors).Execute()

            # set the buffer size
            if self.model_part.GetBufferSize() < self.min_buffer_size:
                self.model_part.SetBufferSize(self.min_buffer_size)

            # initialize the adjoint model part using primal results
            primal_model_part = self.model.GetModelPart(self.primal_model_part_name)
            variable_utils = kratos.VariableUtils()
            variable_utils.CopyModelPartNodalVar(kratos.CONDUCTIVITY, primal_model_part, self.model_part, 0)
            variable_utils.CopyModelPartNodalVar(kratos.TEMPERATURE, primal_model_part, self.model_part, 0)
            variable_utils.CopyModelPartNodalVar(kratos.HEAT_FLUX, primal_model_part, self.model_part, 0)
            variable_utils.CopyModelPartNodalVar(kratos.FACE_HEAT_FLUX, primal_model_part, self.model_part, 0)

            self.ImportMaterials()

    def ImportMaterials(self):
        materials_filename = self.settings["material_import_settings"]["materials_filename"].GetString()
        if (materials_filename != ""):
            with open(materials_filename, 'r') as parameter_file:
                materials = kratos.Parameters(parameter_file.read())

            for i in range(materials["properties"].size()):
                model_part = self.model.GetModelPart(materials["properties"][i]["model_part_name"].GetString())
                mat = materials["properties"][i]["Material"]
                var_utils = kratos.VariableUtils()
                for key, value in mat["Variables"].items():
                    var = kratos.KratosGlobals.GetVariable(key)
                    #if not model_part.HasNodalSolutionStepVariable(var):
                    #     raise Exception("Trying to set variable {0} on nodes, but the variable is not in nodal data.".format(var.Name()))
                    if model_part.HasNodalSolutionStepVariable(var):
                        if value.IsDouble():
                            var_utils.SetScalarVar(var, value.GetDouble(), model_part.Nodes)
                        elif value.IsVector():
                            var_utils.SetVectorVar(var, value.GetVector(), model_part.Nodes)
                        else:
                            raise ValueError("Type of value is not available")
                    else:
                        kratos.Logger.PrintWarning("Ignoring variable {0} given by the materials file, since it is not a nodal variable used by this solver.".format(var.Name()))

    def DefineConvectionDiffusionSettings(self,settings):
        convection_diffusion_settings = kratos.ConvectionDiffusionSettings()

        convection_diffusion_settings.SetDiffusionVariable(
            kratos.KratosGlobals.GetVariable(settings["diffusion_variable"].GetString()))
        convection_diffusion_settings.SetUnknownVariable(
            kratos.KratosGlobals.GetVariable(settings["unknown_variable"].GetString()))
        convection_diffusion_settings.SetVolumeSourceVariable(
            kratos.KratosGlobals.GetVariable(settings["volume_source_variable"].GetString()))
        convection_diffusion_settings.SetSurfaceSourceVariable(
            kratos.KratosGlobals.GetVariable(settings["surface_source_variable"].GetString()))

        self.model_part.ProcessInfo.SetValue(kratos.CONVECTION_DIFFUSION_SETTINGS,convection_diffusion_settings)

    def GetComputingModelPart(self):
        return self.model_part

    def Initialize(self):

        if self.settings["response_function_settings"]["response_type"].GetString() == "point_temperature":
            self.response_function = convdiff.LocalTemperatureAverageResponseFunction(self.settings["response_function_settings"]["custom_settings"],self.model_part)
        else:
            raise Exception("invalid response_type: " + self.settings["response_function_settings"]["response_type"].GetString())

        self.sensitivity_builder = kratos.SensitivityBuilder(self.settings["sensitivity_settings"], self.model_part, self.response_function)

        import KratosMultiphysics.python_linear_solver_factory as linear_solver_factory
        self.linear_solver = linear_solver_factory.ConstructSolver(self.settings["linear_solver_settings"])

        self.time_scheme = kratos.ResidualBasedAdjointStaticScheme(self.response_function)

        builder_and_solver = kratos.ResidualBasedBlockBuilderAndSolver(self.linear_solver)

        self.solver = kratos.ResidualBasedLinearStrategy(self.model_part,
                                                         self.time_scheme,
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

    def AdvanceInTime(self, current_time):
        dt = self.ComputeDeltaTime()
        new_time = current_time + dt
        self.model_part.ProcessInfo[kratos.STEP] += 1
        self.model_part.CloneTimeStep(new_time)

        return new_time

    def ComputeDeltaTime(self):
        return self.settings["time_stepping"]["time_step"].GetDouble()

    def _set_elements_and_conditions(self):

        domain_size = self.model_part.ProcessInfo[kratos.DOMAIN_SIZE]
        comm = self.model_part.GetCommunicator().GetDataCommunicator()

        element_name = self.settings["element_replace_settings"]["element_name"].GetString()
        condition_name = self.settings["element_replace_settings"]["condition_name"].GetString()

        num_nodes_elements = 0
        for elem in self.model_part.Elements:
            num_nodes_elements = len(elem.GetNodes())
            break
        num_nodes_elements = comm.MaxAll(num_nodes_elements)

        if element_name == "AdjointDiffusionElement":
            name_string = "{0}{1}D{2}N".format(element_name,domain_size, num_nodes_elements)
            self.settings["element_replace_settings"]["element_name"].SetString(name_string)

        num_nodes_conditions = 0
        for cond in self.model_part.Conditions:
            num_nodes_conditions = len(cond.GetNodes())
            break
        num_nodes_conditions = comm.MaxAll(num_nodes_conditions)

        if condition_name == "AdjointThermalFace":
            name_string = "{0}{1}D{2}N".format(condition_name,domain_size, num_nodes_conditions)
            self.settings["element_replace_settings"]["condition_name"].SetString(name_string)

        ## Call the replace elements and conditions process
        kratos.ReplaceElementsAndConditionsProcess(self.model_part, self.settings["element_replace_settings"]).Execute()
