"""Load simulation results in the initialization step from HDF5.

license: HDF5Application/license.txt
"""


__all__ = ["Factory"]


import KratosMultiphysics
from KratosMultiphysics.HDF5Application.core.controllers import DefaultController
from KratosMultiphysics.HDF5Application.core.operations.aggregated_operations import ControlledOperation
from KratosMultiphysics.HDF5Application.core.operations.aggregated_operations import AggregatedControlledOperations
from KratosMultiphysics.HDF5Application.core.processes import HDF5Process
from KratosMultiphysics.HDF5Application.core.operations.model_part import *


def Factory(settings: KratosMultiphysics.Parameters, model: KratosMultiphysics.Model) -> HDF5Process:
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise RuntimeError(f"settings should be of type KratosMultiphysics.Parameters")
    if not isinstance(model, KratosMultiphysics.Model):
        raise RuntimeError(f"model should be of type KratosMultiphysics.Model")

    return ImportModelPartFromHDF5Process(model, settings["Parameters"])


class ImportModelPartFromHDF5Process(HDF5Process):
    @classmethod
    def GetDefaultParameters(cls) -> KratosMultiphysics.Parameters:
        return KratosMultiphysics.Parameters("""
            {
                "model_part_name": "PLEASE_SPECIFY_MODEL_PART_NAME",
                "file_settings": {
                    "file_name"        : "<model_part_name>",
                    "time_format"      : "0.4f",
                    "file_access_mode" : "read_only",
                    "echo_level"       :  0
                },
                "model_part_input_settings": {
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
                "condition_flag_value_settings": {
                    "prefix"           : "/ResultsData/ConditionFlagValues/",
                    "list_of_variables": [],
                    "time_format"      : "0.4f"
                }
            }""")

    def __init__(self, model: KratosMultiphysics.Model, parameters: KratosMultiphysics.Parameters) -> None:
        parameters.RecursivelyValidateAndAssignDefaults(self.GetDefaultParameters())
        model_part = model[parameters["model_part_name"].GetString()]
        super().__init__()

        # create default controller
        default_controller = DefaultController()

        # create the aggregated operation with hdf5 file settings
        operations = AggregatedControlledOperations(model_part, parameters["file_settings"])

        # now adding temporal outputs.
        operations.AddControlledOperation(ControlledOperation(ModelPartInput, self._GetValidatedParameters("model_part_input_settings", parameters), default_controller))
        operations.AddControlledOperation(ControlledOperation(NodalSolutionStepDataInput, self._GetValidatedParameters("nodal_solution_step_data_settings", parameters), default_controller))
        operations.AddControlledOperation(ControlledOperation(NodalDataValueInput, self._GetValidatedParameters("nodal_data_value_settings", parameters), default_controller))
        operations.AddControlledOperation(ControlledOperation(NodalFlagValueInput, self._GetValidatedParameters("nodal_flag_value_settings", parameters), default_controller))
        operations.AddControlledOperation(ControlledOperation(ElementDataValueInput, self._GetValidatedParameters("element_data_value_settings", parameters), default_controller))
        operations.AddControlledOperation(ControlledOperation(ElementFlagValueInput, self._GetValidatedParameters("element_flag_value_settings", parameters), default_controller))
        operations.AddControlledOperation(ControlledOperation(ConditionDataValueInput, self._GetValidatedParameters("condition_data_value_settings", parameters), default_controller))
        operations.AddControlledOperation(ControlledOperation(ConditionFlagValueInput, self._GetValidatedParameters("condition_flag_value_settings", parameters), default_controller))

        # now add all operations to PrintOutput method
        self.AddInitialize(operations)