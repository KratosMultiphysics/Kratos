from shutil import copytree
from pathlib import Path
from importlib import import_module

import KratosMultiphysics as Kratos
import KratosMultiphysics.HDF5Application as KratosHDF5
from KratosMultiphysics.analysis_stage import AnalysisStage
from KratosMultiphysics.project import Project
from KratosMultiphysics.orchestrators.orchestrator import Orchestrator
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView
from KratosMultiphysics.OptimizationApplication.execution_policies.execution_policy import ExecutionPolicy
from KratosMultiphysics.OptimizationApplication.utilities.helper_utilities import GetClassModuleFromKratos, WorkFolderScope
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem_utilities import GetComponentHavingDataByFullName, GetComponentValueByFullName
from KratosMultiphysics.HDF5Application.core.file_io import OpenHDF5File

def Factory(model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem) -> ExecutionPolicy:
    if not parameters.Has("name"):
        raise RuntimeError(f"IndependentAnalysisExecutionPolicy instantiation requires a \"name\" in parameters [ parameters = {parameters}].")
    if not parameters.Has("settings"):
        raise RuntimeError(f"IndependentAnalysisExecutionPolicy instantiation requires a \"settings\" in parameters [ parameters = {parameters}].")

    return IndependentAnalysisExecutionPolicy(parameters["name"].GetString(), model, parameters["settings"], optimization_problem)

class IndependentAnalysisExecutionPolicy(ExecutionPolicy):
    def __init__(self, name: str, model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
        super().__init__(name)

        self.model = model
        self.optimization_problem = optimization_problem

        default_settings = Kratos.Parameters("""{
            "analysis_module"         : "KratosMultiphysics",
            "analysis_type"           : "",
            "analysis_model_part_name": "",
            "analysis_settings"       : {
                "input_folder"                : "",
                "output_folder"               : "Optimization_Results/<name>/<step>",
                "project_parameters_file_name": "",
                "output_settings"             : {
                    "hdf5_file_name": "",
                    "list_of_values": [
                        {
                            "component_name": "",
                            "value_name"    : ""
                        }
                    ]
                }
            }
        }""")

        parameters.RecursivelyValidateAndAssignDefaults(default_settings)

        self.analysis_module = parameters["analysis_module"].GetString()
        self.analysis_type = parameters["analysis_type"].GetString()
        self.analysis_model_part_name = parameters["analysis_model_part_name"].GetString()

        analysis_settings = parameters["analysis_settings"]
        self.analysis_input_folder = analysis_settings["input_folder"].GetString()
        self.analysis_output_folder = analysis_settings["output_folder"].GetString()
        self.analysis_project_parameters_file_name = analysis_settings["project_parameters_file_name"].GetString()

        self.analysis_output_h5_file_name = analysis_settings["output_settings"]["hdf5_file_name"].GetString()
        self.analysis_list_of_values: 'list[tuple[str, str]]' = []
        for value_settings in analysis_settings["output_settings"]["list_of_values"].values():
            value_settings.ValidateAndAssignDefaults(default_settings["analysis_settings"]["output_settings"]["list_of_values"].values()[0])
            self.analysis_list_of_values.append((value_settings["component_name"].GetString(), value_settings["value_name"].GetString()))

        if self.analysis_module == "KratosMultiphysics":
            self.analysis_full_module, self.analysis_type = GetClassModuleFromKratos(self.analysis_type)
        else:
            self.analysis_full_module = f"{self.analysis_module}.{Kratos.StringUtilities.ConvertCamelCaseToSnakeCase(self.analysis_type)}"

        self.current_analysis = None

    def Initialize(self) -> None:
        pass

    def Check(self) -> None:
        pass

    def Finalize(self) -> None:
        pass

    def Execute(self):
        self.analysis_path = self.analysis_output_folder.replace("<name>", self.GetName())
        self.analysis_path = self.analysis_path.replace("<step>", str(self.optimization_problem.GetStep()))
        self.analysis_path = Path(self.analysis_path)
        self.analysis_path.mkdir(parents=True, exist_ok=True)

        copytree(self.analysis_input_folder, self.analysis_path, dirs_exist_ok=True)

        with WorkFolderScope(self.analysis_path):
            # create a model always.
            analysis_model = Kratos.Model()

            hdf5_file_parameter_settings = Kratos.Parameters("""{
                "file_name"       : """ + f"\"{self.analysis_output_h5_file_name}\"" + """,
                "file_access_mode": "truncate"
            }""")
            with OpenHDF5File(hdf5_file_parameter_settings, self.model[self.analysis_model_part_name]) as h5_file:
                attributes = Kratos.Parameters("""{}""")

                for component_name, value_name in self.analysis_list_of_values:
                    component = GetComponentHavingDataByFullName(component_name, self.optimization_problem)
                    value = GetComponentValueByFullName(ComponentDataView(component, self.optimization_problem), value_name)

                    h5_path = f"/{component_name.replace(".", "/")}"
                    h5_file.AddPath(h5_path)
                    h5_value_name = value_name.split(":")[0]
                    if isinstance(value, int):
                        attributes.AddInt(h5_value_name, value)
                    elif isinstance(value, float):
                        attributes.AddDouble(h5_value_name, value)
                    elif isinstance(value, bool):
                        attributes.AddBool(h5_value_name, value)
                    elif isinstance(value, (Kratos.TensorAdaptors.DoubleTensorAdaptor, Kratos.TensorAdaptors.IntTensorAdaptor, Kratos.TensorAdaptors.BoolTensorAdaptor)):
                        ta_parameters = Kratos.Parameters("""{
                            "prefix": """ + f"\"{h5_path}/\"" + """
                        }""")
                        KratosHDF5.TensorAdaptorIO(ta_parameters, h5_file).Write(h5_value_name, value)
                    else:
                        raise RuntimeError(f"Unsupported value type [ component_name: {component_name}, value_name: {value_name}, value = {value} ].")
                if not attributes.IsEquivalentTo(Kratos.Parameters("""{}""")):
                    h5_file.WriteAttribute(h5_path, attributes)

            analysis_type = getattr(import_module(self.analysis_full_module), self.analysis_type)
            with open(self.analysis_project_parameters_file_name, "r") as file_input:
                analysis_parameters = Kratos.Parameters(file_input.read())
            if AnalysisStage in analysis_type.mro():
                # the analysis type is derived from AnalysisStage
                self.current_analysis: AnalysisStage = getattr(import_module(self.analysis_full_module), self.analysis_type)(analysis_model, analysis_parameters)
            elif Orchestrator in analysis_type.mro():
                # the analysis type is derive from the Orchestrator
                project = Project(analysis_parameters)
                self.current_analysis: Orchestrator = getattr(import_module(self.analysis_full_module), self.analysis_type)(project)

            self.current_analysis.Run()

    def GetPath(self) -> str:
        return self.analysis_path

    def GetAnalysisModelPart(self):
        if self.current_analysis is not None:
            if isinstance(self.current_analysis, AnalysisStage):
                return self.current_analysis._GetSolver().GetComputingModelPart()
            elif isinstance(self.current_analysis, Orchestrator):
                return self.current_analysis.GetProject().GetModel()[self.analysis_model_part_name]
            else:
                raise RuntimeError(f"Unsupported analysis type = {self.current_analysis}.")
        else:
            return self.model[self.analysis_model_part_name]
