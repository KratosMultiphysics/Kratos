"""Store temporal simulation results for a single mesh with HDF5.

This process:
 - stores the initial model part in an .h5 file.
 - stores historical and non-historical results in one .h5 file per output step.

This process works with or without MPI.

license: HDF5Application/license.txt
"""


__all__ = ["Factory"]


import KratosMultiphysics
from KratosMultiphysics.HDF5Application.core.controllers import SingleTimeController
from KratosMultiphysics.HDF5Application.core.operations.aggregated_operations import ControlledOperation
from KratosMultiphysics.HDF5Application.core.operations.aggregated_operations import AggregatedControlledOperations
from KratosMultiphysics.HDF5Application.core.processes import HDF5OutputProcess
from KratosMultiphysics.HDF5Application.core.operations.model_part import *


def Factory(settings: KratosMultiphysics.Parameters, model: KratosMultiphysics.Model) -> HDF5OutputProcess:
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise RuntimeError(f"settings should be of type KratosMultiphysics.Parameters")
    if not isinstance(model, KratosMultiphysics.Model):
        raise RuntimeError(f"model should be of type KratosMultiphysics.Model")

    return SingleMeshTemporalOutputProcess(model, settings["Parameters"])


class SingleMeshTemporalOutputProcess(HDF5OutputProcess):
    @classmethod
    def GetDefaultParameters(cls) -> KratosMultiphysics.Parameters:
        return KratosMultiphysics.Parameters("""
            {
                "model_part_name": "PLEASE_SPECIFY_MODEL_PART_NAME",
                "file_settings": {
                    "file_name"        : "<model_part_name>-<time>.h5",
                    "time_format"      : "0.4f",
                    "file_access_mode" : "exclusive",
                    "max_files_to_keep": "unlimited",
                    "echo_level"       :  0
                },
                "output_time_settings": {
                    "output_control_type": "time",
                    "output_interval"    : 1.0
                },
                "model_part_output_settings": {
                    "prefix"     : "/ModelData",
                    "time_format": "0.4f"
                },
                "nodal_solution_step_data_settings": {
                    "prefix"           : "/ResultsData/NodalSolutionStepData/",
                    "list_of_variables": [],
                    "time_format"      : "0.4f"
                },
                "nodal_data_value_settings": {
                    "prefix"           : "/ResultsData/NodalDataValues/",
                    "list_of_variables": [],
                    "time_format"      : "0.4f"
                },
                "nodal_flag_value_settings": {
                    "prefix"           : "/ResultsData/NodalFlagValues/",
                    "list_of_variables": [],
                    "time_format"      : "0.4f"
                },
                "element_data_value_settings": {
                    "prefix"           : "/ResultsData/ElementDataValues/",
                    "list_of_variables": [],
                    "time_format"      : "0.4f"
                },
                "element_gauss_point_value_settings": {
                    "prefix"           : "/ResultsData/ElementGaussPointValues/",
                    "list_of_variables": [],
                    "time_format"      : "0.4f"
                },
                "element_flag_value_settings": {
                    "prefix"           : "/ResultsData/ElementFlagValues/",
                    "list_of_variables": [],
                    "time_format"      : "0.4f"
                },
                "condition_data_value_settings": {
                    "prefix"           : "/ResultsData/ConditionDataValues/",
                    "list_of_variables": [],
                    "time_format"      : "0.4f"
                },
                "condition_gauss_point_value_settings": {
                    "prefix"           : "/ResultsData/ConditionGaussPointValues/",
                    "list_of_variables": [],
                    "time_format"      : "0.4f"
                },
                "condition_flag_value_settings": {
                    "prefix"           : "/ResultsData/ConditionFlagValues/",
                    "list_of_variables": [],
                    "time_format"      : "0.4f"
                }
            }""")

    def __init__(self, model: KratosMultiphysics.Model, parameters: KratosMultiphysics.Parameters) -> None:
        super().__init__()
        parameters.RecursivelyValidateAndAssignDefaults(self.GetDefaultParameters())

        model_part = model[parameters["model_part_name"].GetString()]

        # create temporal controller
        temporal_controller_settings = parameters["output_time_settings"]
        temporal_controller_settings.AddString("model_part_name", model_part.FullName())
        temporal_controller = KratosMultiphysics.TemporalController(model, temporal_controller_settings)

        # create the aggregated operation with hdf5 file settings
        operations = AggregatedControlledOperations(model_part, parameters["file_settings"])

        # adding one time mesh output
        operations.AddControlledOperation(ControlledOperation(ModelPartOutput, parameters["model_part_output_settings"], SingleTimeController(temporal_controller)))

        # now adding temporal outputs.
        operations.AddControlledOperation(ControlledOperation(NodalSolutionStepDataOutput, parameters["nodal_solution_step_data_settings"], temporal_controller))
        operations.AddControlledOperation(ControlledOperation(NodalDataValueOutput, parameters["nodal_data_value_settings"], temporal_controller))
        operations.AddControlledOperation(ControlledOperation(NodalFlagValueOutput, parameters["nodal_flag_value_settings"], temporal_controller))
        operations.AddControlledOperation(ControlledOperation(ElementDataValueOutput, parameters["element_data_value_settings"], temporal_controller))
        operations.AddControlledOperation(ControlledOperation(ElementGaussPointOutput, parameters["element_gauss_point_value_settings"], temporal_controller))
        operations.AddControlledOperation(ControlledOperation(ElementFlagValueOutput, parameters["element_flag_value_settings"], temporal_controller))
        operations.AddControlledOperation(ControlledOperation(ConditionDataValueOutput, parameters["condition_data_value_settings"], temporal_controller))
        operations.AddControlledOperation(ControlledOperation(ConditionGaussPointOutput, parameters["condition_gauss_point_value_settings"], temporal_controller))
        operations.AddControlledOperation(ControlledOperation(ConditionFlagValueOutput, parameters["condition_flag_value_settings"], temporal_controller))

        # now add all operations to PrintOutput method
        self.AddPrintOutput(operations)