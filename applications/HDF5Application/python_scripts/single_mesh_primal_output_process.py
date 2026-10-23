"""Store primal simulation results for a single mesh with HDF5.

This process:
 - stores the initial model part in an .h5 file.
 - stores historical and non-historical results in one .h5 file per output step.
   The variable ACCELERATION is stored in Bossak weighted form for
   time-dependent adjoint simulations.

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

    return SingleMeshPrimalOutputProcess(model, settings["Parameters"])


class SingleMeshPrimalOutputProcess(HDF5OutputProcess):
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
                    "time_format": "0.4f",
                    "custom_attributes": {}
                },
                "nodal_solution_step_data_settings": {
                    "prefix"           : "/ResultsData/NodalSolutionStepData/",
                    "list_of_variables": [],
                    "time_format"      : "0.4f",
                    "custom_attributes": {},
                    "alpha_bossak"     : -0.3
                },
                "nodal_data_value_settings": {
                    "prefix"           : "/ResultsData/NodalDataValues/",
                    "list_of_variables": [],
                    "time_format"      : "0.4f",
                    "custom_attributes": {}
                },
                "nodal_flag_value_settings": {
                    "prefix"           : "/ResultsData/NodalFlagValues/",
                    "list_of_variables": [],
                    "time_format"      : "0.4f",
                    "custom_attributes": {}
                },
                "element_data_value_settings": {
                    "prefix"           : "/ResultsData/ElementDataValues/",
                    "list_of_variables": [],
                    "time_format"      : "0.4f",
                    "custom_attributes": {}
                },
                "element_gauss_point_value_settings": {
                    "prefix"           : "/ResultsData/ElementGaussPointValues/",
                    "list_of_variables": [],
                    "time_format"      : "0.4f",
                    "custom_attributes": {}
                },
                "element_flag_value_settings": {
                    "prefix"           : "/ResultsData/ElementFlagValues/",
                    "list_of_variables": [],
                    "time_format"      : "0.4f",
                    "custom_attributes": {}
                },
                "condition_data_value_settings": {
                    "prefix"           : "/ResultsData/ConditionDataValues/",
                    "list_of_variables": [],
                    "time_format"      : "0.4f",
                    "custom_attributes": {}
                },
                "condition_gauss_point_value_settings": {
                    "prefix"           : "/ResultsData/ConditionGaussPointValues/",
                    "list_of_variables": [],
                    "time_format"      : "0.4f",
                    "custom_attributes": {}
                },
                "condition_flag_value_settings": {
                    "prefix"           : "/ResultsData/ConditionFlagValues/",
                    "list_of_variables": [],
                    "time_format"      : "0.4f",
                    "custom_attributes": {}
                }
            }""")

    def __init__(self, model: KratosMultiphysics.Model, parameters: KratosMultiphysics.Parameters) -> None:
        parameters.ValidateAndAssignDefaults(self.GetDefaultParameters())
        model_part = model[parameters["model_part_name"].GetString()]
        super().__init__()

        # create temporal controller
        temporal_controller_settings = parameters["output_time_settings"]
        temporal_controller_settings.AddString("model_part_name", model_part.FullName())
        temporal_controller = KratosMultiphysics.OutputController(model, temporal_controller_settings)

        # create the aggregated operation with hdf5 file settings
        operations = AggregatedControlledOperations(model_part, self._GetValidatedParameters("file_settings", parameters))

        # adding one time mesh output
        operations.AddControlledOperation(ControlledOperation(ModelPartOutput, self._GetOperationParameters("model_part_output_settings", parameters), SingleTimeController(temporal_controller)))

        # now adding temporal outputs.
        operations.AddControlledOperation(ControlledOperation(PrimalBossakOutput, self._GetOperationParameters("nodal_solution_step_data_settings", parameters), temporal_controller))
        operations.AddControlledOperation(ControlledOperation(NodalDataValueOutput, self._GetOperationParameters("nodal_data_value_settings", parameters), temporal_controller))
        operations.AddControlledOperation(ControlledOperation(NodalFlagValueOutput, self._GetOperationParameters("nodal_flag_value_settings", parameters), temporal_controller))
        operations.AddControlledOperation(ControlledOperation(ElementDataValueOutput, self._GetOperationParameters("element_data_value_settings", parameters), temporal_controller))
        operations.AddControlledOperation(ControlledOperation(ElementGaussPointOutput, self._GetOperationParameters("element_gauss_point_value_settings", parameters), temporal_controller))
        operations.AddControlledOperation(ControlledOperation(ElementFlagValueOutput, self._GetOperationParameters("element_flag_value_settings", parameters), temporal_controller))
        operations.AddControlledOperation(ControlledOperation(ConditionDataValueOutput, self._GetOperationParameters("condition_data_value_settings", parameters), temporal_controller))
        operations.AddControlledOperation(ControlledOperation(ConditionGaussPointOutput, self._GetOperationParameters("condition_gauss_point_value_settings", parameters), temporal_controller))
        operations.AddControlledOperation(ControlledOperation(ConditionFlagValueOutput, self._GetOperationParameters("condition_flag_value_settings", parameters), temporal_controller))

        # now add all operations to PrintOutput method
        self.AddPrintOutput(operations)
