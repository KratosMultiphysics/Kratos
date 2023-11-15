from pathlib import Path

import KratosMultiphysics as Kratos
import KratosMultiphysics.HDF5Application as KratosHDF5
import KratosMultiphysics.OptimizationApplication as KratosOA

from KratosMultiphysics.HDF5Application.core.file_io import OpenHDF5File
from KratosMultiphysics.DigitalTwinApplication.utilities.data_utils import SupportedVariableUnionType
from KratosMultiphysics.DigitalTwinApplication.utilities.expression_utils import GetContainerExpression
from KratosMultiphysics.DigitalTwinApplication.utilities.expression_utils import GetContainerExpressionType
from KratosMultiphysics.DigitalTwinApplication.utilities.expression_utils import ExpressionDataLocation
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem

class DataIO:
    def __init__(self, output_path: str, model_part: Kratos.ModelPart, optimization_problem: OptimizationProblem) -> None:
        self.model_part = model_part
        self.output_path = Path(output_path)
        self.list_of_data_ios: 'list[tuple[ExpressionDataLocation, SupportedVariableUnionType]]' = []
        self.optimization_problem = optimization_problem

    def AddContainerVariable(self, container_type: ExpressionDataLocation, variable: SupportedVariableUnionType) -> None:
        self.list_of_data_ios.append((container_type, variable))

    def Read(self) -> None:
        if len(self.list_of_data_ios) == 0:
            return

        h5_file_params = Kratos.Parameters("""{
            "file_name": "",
            "file_access_mode": "read_only"
        }""")
        h5_file_params["file_name"].SetString(str(self.output_path / f"{self.model_part.FullName()}_{self.optimization_problem.GetStep()}.h5"))
        with OpenHDF5File(h5_file_params, self.model_part) as h5_file:
            h5_data_params = Kratos.Parameters("""{
                "prefix": "/SystemIdentificationData"
            }""")
            for container_type, variable in self.list_of_data_ios:
                h5_data_params["prefix"].SetString(f"/SystemIdentificationData/{container_type.name}/")
                h5_io = KratosHDF5.ExpressionIO(h5_data_params, h5_file)
                exp = GetContainerExpressionType(container_type)(self.model_part)
                h5_io.Read(variable.Name(), exp)

                if container_type == ExpressionDataLocation.NodeHistorical:
                    Kratos.Expression.VariableExpressionIO.Write(exp, variable, True)
                elif container_type == ExpressionDataLocation.NodeNonHistorical:
                    Kratos.Expression.VariableExpressionIO.Write(exp, variable, False)
                elif container_type in [ExpressionDataLocation.Condition, ExpressionDataLocation.Element]:
                    Kratos.Expression.VariableExpressionIO.Write(exp, variable)
                elif container_type in [ExpressionDataLocation.ConditionProperties, ExpressionDataLocation.ElementProperties]:
                    KratosOA.PropertiesVariableExpressionIO.Write(exp, variable)
                else:
                    raise RuntimeError(f"Unsupported container type = {container_type.name}.")

    def Write(self) -> None:
        if len(self.list_of_data_ios) == 0:
            return

        h5_file_params = Kratos.Parameters("""{
            "file_name": "",
            "file_access_mode": "truncate"
        }""")
        h5_file_params["file_name"].SetString(str(self.output_path / f"{self.model_part.FullName()}_{self.optimization_problem.GetStep()}.h5"))

        with OpenHDF5File(h5_file_params, self.model_part) as h5_file:
            h5_data_params = Kratos.Parameters("""{
                "prefix": "/SystemIdentificationData"
            }""")
            for container_type, variable in self.list_of_data_ios:
                h5_data_params["prefix"].SetString(f"/SystemIdentificationData/{container_type.name}/")
                h5_io = KratosHDF5.ExpressionIO(h5_data_params, h5_file)
                exp = GetContainerExpression(self.model_part, container_type, variable)
                h5_io.Write(variable.Name(), exp)


def GetDataIO(model: Kratos.Model, output_path: str, params: Kratos.Parameters, optimization_problem: OptimizationProblem) -> DataIO:
    default_params = Kratos.Parameters("""{
        "model_part_name"              : "PLEASE_SPECIFY_MODEL_PART_NAME",
        "nodal_hist_variable_names"    : [],
        "nodal_non_hist_variable_names": [],
        "condition_variable_names"     : [],
        "element_variable_names"       : [],
        "condition_property_names"     : [],
        "element_property_names"       : []
    }""")
    params.ValidateAndAssignDefaults(default_params)

    nodal_hist_variables = [Kratos.KratosGlobals.GetVariable(var_name) for var_name in params["nodal_hist_variable_names"].GetStringArray()]
    nodal_non_hist_variables = [Kratos.KratosGlobals.GetVariable(var_name) for var_name in params["nodal_non_hist_variable_names"].GetStringArray()]
    condition_variables = [Kratos.KratosGlobals.GetVariable(var_name) for var_name in params["condition_variable_names"].GetStringArray()]
    element_variables = [Kratos.KratosGlobals.GetVariable(var_name) for var_name in params["element_variable_names"].GetStringArray()]
    condition_property_variables = [Kratos.KratosGlobals.GetVariable(var_name) for var_name in params["condition_property_names"].GetStringArray()]
    element_property_variables = [Kratos.KratosGlobals.GetVariable(var_name) for var_name in params["element_property_names"].GetStringArray()]

    model_part = model[params["model_part_name"].GetString()]

    data_io = DataIO(output_path, model_part, optimization_problem)
    for var in nodal_hist_variables:
        data_io.AddContainerVariable(ExpressionDataLocation.NodeHistorical, var)
    for var in nodal_non_hist_variables:
        data_io.AddContainerVariable(ExpressionDataLocation.NodeNonHistorical, var)
    for var in condition_variables:
        data_io.AddContainerVariable(ExpressionDataLocation.Condition, var)
    for var in element_variables:
        data_io.AddContainerVariable(ExpressionDataLocation.Element, var)
    for var in condition_property_variables:
        data_io.AddContainerVariable(ExpressionDataLocation.ConditionProperties, var)
    for var in element_property_variables:
        data_io.AddContainerVariable(ExpressionDataLocation.ElementProperties, var)

    return data_io