from pathlib import Path
import typing
import KratosMultiphysics as Kratos
import KratosMultiphysics.SystemIdentificationApplication as KratosSI
import KratosMultiphysics.HDF5Application as KratosHDF5
import KratosMultiphysics.OptimizationApplication as KratosOA

from KratosMultiphysics.HDF5Application.core.file_io import OpenHDF5File
from KratosMultiphysics.HDF5Application.single_mesh_temporal_output_process import SingleMeshTemporalOutputProcessFactory

from KratosMultiphysics.SystemIdentificationApplication.utilities.sensor_utils import CreateSensors, GetSensors
from KratosMultiphysics.SystemIdentificationApplication.utilities.sensor_utils import PrintSensorListToCSV
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView
from KratosMultiphysics.OptimizationApplication.responses.response_function import ResponseFunction


def Factory(arg1: Kratos.Parameters, arg2: Kratos.Model, optimization_problem: 'typing.Optional[OptimizationProblem]' = None):
    if isinstance(arg1, Kratos.Parameters) and isinstance(arg2, Kratos.Model):
        # this is to use PrimalHDF5OutputProcess as a normal kratos process
        return PrimalHDF5OutputProcess(arg2, arg1["Parameters"], optimization_problem)
    elif isinstance(arg1, Kratos.Model) and isinstance(arg2, Kratos.Parameters):
        # this is to use PrimalHDF5OutputProcess as a process for optimization application where optimization_problem cannot be none
        return PrimalHDF5OutputProcess(arg1, arg2["settings"], optimization_problem)
    else:
        raise RuntimeError("Argument mismatch.")

class PrimalHDF5OutputProcess(Kratos.OutputProcess):
    def __init__(self, model: Kratos.Model, settings: Kratos.Parameters, optimization_problem: 'typing.Optional[OptimizationProblem]' = None):
        Kratos.OutputProcess.__init__(self)

        default_settings = Kratos.Parameters("""{
            "model_part_name": "PLEASE_PROVIDE_A_MODEL_PART_NAME",
            "nodal_solution_step_variables": ["DISPLACEMENT", "ACCELERATION"],
            "output_file_name": "hdf5_output/<model_part_name>_T_<step>.h5"
        }""")

        if optimization_problem is not None:
            default_settings.AddEmptyValue("response_name")
            default_settings["response_name"].SetString("")
        settings.ValidateAndAssignDefaults(default_settings)

        self.model_part = model[settings["model_part_name"].GetString()]
        self.output_file_name = Path(settings["output_file_name"].GetString())
        self.optimization_problem = optimization_problem

        self.response_name = self.model_part.GetValue(KratosOA.EXECUTION_POLICY_NAME)
        self.nodal_solution_step_variables = settings["nodal_solution_step_variables"].GetStringArray()


    def PrintOutput(self) -> None:
        s_name = str(self.output_file_name)
        s_name = s_name.replace("<model_part_name>", f"{self.model_part.Name}")
        s_name_orig = s_name
        if self.optimization_problem is None:
            s_name = s_name.replace("<step>", f"{self.model_part.ProcessInfo[Kratos.STEP]}")
        else:
            s_name = s_name.replace("<step>", f"{self.optimization_problem.GetStep()}")
       
        #print("s_name is ", s_name)
        params = Kratos.Parameters("""{
            "file_access_mode": "truncate"
        }""")
        params.AddString("file_name", s_name)

        prefix_settings = Kratos.Parameters("""{                         
        }""")
        if not prefix_settings.Has("list_of_variables"):
            prefix_settings.AddStringArray("list_of_variables", self.nodal_solution_step_variables)
        if not prefix_settings.Has("prefix"):
            prefix = "/ResultsData"
            prefix_settings.AddString("prefix", prefix)

        with OpenHDF5File(params, self.model_part) as h5_file:
            expio = KratosHDF5.HDF5NodalSolutionStepDataIO(prefix_settings, h5_file)
            expio.Write(self.model_part)

        exec_policy_name = self.model_part.GetValue(KratosOA.EXECUTION_POLICY_NAME)

        KratosOA.OptAppModelPartUtils.LogModelPartStatus(self.model_part, f"{exec_policy_name}_PRIMAL_DATA_written")
        KratosOA.OptAppModelPartUtils.LogModelPartStatus(self.model_part, f"{exec_policy_name}_PRIMAL_DATA_Path:{s_name_orig}:{prefix}:{self.nodal_solution_step_variables}")


        




