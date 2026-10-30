"""In-process execution backend: runs the AnalysisStage in this interpreter.

Convenience mode for small problems: no serialization overhead, but the solve
shares the process (and the GIL) with the training loop. For MPI-parallel
Kratos runs or HPC scheduling use the subprocess backend instead.
"""

import os
import string
from importlib import import_module
from pathlib import Path

import numpy

import KratosMultiphysics as Kratos
from KratosMultiphysics.PhysicsNeMoApplication.active_learning.execution_backends.base_backend import KratosExecutionBackend
from KratosMultiphysics.PhysicsNeMoApplication.active_learning.sample_io import KratosALSample, ApplyParameterOverrides
from KratosMultiphysics.PhysicsNeMoApplication.utilities.tensor_adaptor_dataset_utils import GetTensorAdaptor


class InProcessBackend(KratosExecutionBackend):
    """Runs an AnalysisStage per sample inside the current process."""

    def __init__(self, settings: Kratos.Parameters) -> None:
        default_settings = Kratos.Parameters("""{
            "project_parameters_file" : "ProjectParameters.json",
            "analysis_stage_module"   : "PLEASE_SPECIFY_ANALYSIS_STAGE_MODULE",
            "analysis_stage_name"     : "",
            "working_directory"       : "physics_nemo_al_cases",
            "model_part_name"         : "PLEASE_SPECIFY_MODEL_PART_NAME",
            "output_field_specs"      : [
                {
                    "variable_name" : "PLEASE_SPECIFY_VARIABLE_NAME",
                    "data_location" : "node_historical"
                }
            ]
        }""")
        settings.ValidateAndAssignDefaults(default_settings)
        for i in range(settings["output_field_specs"].size()):
            settings["output_field_specs"][i].ValidateAndAssignDefaults(default_settings["output_field_specs"][0])

        with open(settings["project_parameters_file"].GetString(), "r") as f:
            self.base_parameters = Kratos.Parameters(f.read())
        self.analysis_stage_module = settings["analysis_stage_module"].GetString()
        self.analysis_stage_name = settings["analysis_stage_name"].GetString()
        self.working_directory = Path(settings["working_directory"].GetString())
        self.model_part_name = settings["model_part_name"].GetString()
        self.output_field_specs = [
            (settings["output_field_specs"][i]["variable_name"].GetString(),
             settings["output_field_specs"][i]["data_location"].GetString())
            for i in range(settings["output_field_specs"].size())
        ]

    @property
    def is_external(self) -> bool:
        return False

    def RunCase(self, sample: KratosALSample) -> KratosALSample:
        case_parameters = self.base_parameters.Clone()
        ApplyParameterOverrides(case_parameters, sample.parameters)

        case_directory = self.working_directory / sample.sample_id
        case_directory.mkdir(parents=True, exist_ok=True)

        model = Kratos.Model()
        analysis = _CreateAnalysisStage(self.analysis_stage_module, self.analysis_stage_name, model, case_parameters)

        # Kratos IO is CWD-relative: isolate each sample in its own directory.
        previous_directory = os.getcwd()
        os.chdir(case_directory)
        try:
            analysis.Run()
        finally:
            os.chdir(previous_directory)

        model_part = model[self.model_part_name]
        for variable_name, data_location in self.output_field_specs:
            variable = Kratos.KratosGlobals.GetVariable(variable_name)
            tensor_adaptor = GetTensorAdaptor(model_part, data_location, variable)
            sample.fields[f"{variable_name}__{data_location}"] = numpy.array(tensor_adaptor.data)
        return sample


def _CreateAnalysisStage(module_name: str, stage_name: str, model: Kratos.Model, parameters: Kratos.Parameters):
    """Resolves and instantiates the analysis stage from its module name.

    Follows the CoSimulationApplication kratos_base_wrapper idiom: a module
    "Create" factory wins; otherwise the class named after the module in
    PascalCase; otherwise the explicitly given stage_name.
    """
    analysis_stage_module = import_module(module_name)
    if hasattr(analysis_stage_module, "Create"):
        return analysis_stage_module.Create(model, parameters)

    if not stage_name:
        file_name = module_name.split(".")[-1]
        stage_name = string.capwords(file_name.replace("_", " ")).replace(" ", "")
    if not hasattr(analysis_stage_module, stage_name):
        raise AttributeError(
            f"Module \"{module_name}\" has neither a \"Create\" factory nor a \"{stage_name}\" class; "
            "provide \"analysis_stage_name\" explicitly.")
    return getattr(analysis_stage_module, stage_name)(model, parameters)
