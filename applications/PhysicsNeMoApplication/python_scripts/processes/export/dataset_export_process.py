"""Process exporting model part fields as numpy datasets for ML training.

Writes one ``.npz`` file per sampled step, containing one array per requested
(variable, location) pair plus TIME/STEP metadata. The output is plain numpy —
this module never imports torch or physicsnemo — and is intended as a simple
on-disk format consumable by ``numpy.load`` and, downstream, by
``physicsnemo.datapipes``-style dataset readers.
"""

from pathlib import Path

import numpy

import KratosMultiphysics as Kratos
from KratosMultiphysics.PhysicsNeMoApplication.utilities.tensor_adaptor_dataset_utils import GetTensorAdaptor


def Factory(settings: Kratos.Parameters, model: Kratos.Model) -> "DatasetExportProcess":
    if not isinstance(settings, Kratos.Parameters):
        raise TypeError("Expected input shall be a Parameters object, encapsulating a json string")
    return DatasetExportProcess(model, settings["Parameters"])


class DatasetExportProcess(Kratos.Process):
    """Exports the configured fields every ``output_interval`` steps."""

    def __init__(self, model: Kratos.Model, settings: Kratos.Parameters) -> None:
        super().__init__()

        default_settings = Kratos.Parameters("""{
            "model_part_name" : "PLEASE_SPECIFY_MODEL_PART_NAME",
            "list_of_fields"  : [
                {
                    "variable_name" : "PLEASE_SPECIFY_VARIABLE_NAME",
                    "data_location" : "node_historical"
                }
            ],
            "output_path"     : "physics_nemo_dataset",
            "file_prefix"     : "sample",
            "output_interval" : 1
        }""")
        settings.ValidateAndAssignDefaults(default_settings)
        for i in range(settings["list_of_fields"].size()):
            settings["list_of_fields"][i].ValidateAndAssignDefaults(default_settings["list_of_fields"][0])

        self.model_part = model[settings["model_part_name"].GetString()]
        self.field_specs = [
            (settings["list_of_fields"][i]["variable_name"].GetString(),
             settings["list_of_fields"][i]["data_location"].GetString())
            for i in range(settings["list_of_fields"].size())
        ]
        self.output_path = Path(settings["output_path"].GetString())
        self.file_prefix = settings["file_prefix"].GetString()
        self.output_interval = settings["output_interval"].GetInt()

        if self.output_interval < 1:
            raise ValueError(f"\"output_interval\" must be >= 1 [ output_interval = {self.output_interval} ].")

    def ExecuteInitialize(self) -> None:
        self.output_path.mkdir(parents=True, exist_ok=True)

    def ExecuteFinalizeSolutionStep(self) -> None:
        step = self.model_part.ProcessInfo[Kratos.STEP]
        if step % self.output_interval != 0:
            return

        arrays = {
            "TIME": numpy.array(self.model_part.ProcessInfo[Kratos.TIME]),
            "STEP": numpy.array(step),
        }
        if self.model_part.IsDistributed():
            # Gather ghost-free per-rank blocks to rank 0, which writes ONE
            # file with the entities sorted by global id - identical layout
            # to a serial run, so consumers are rank-count-agnostic.
            from KratosMultiphysics.PhysicsNeMoApplication.distributed import distributed_utils
            is_writing_rank = self.model_part.GetCommunicator().GetDataCommunicator().Rank() == 0
            for variable_name, data_location in self.field_specs:
                _, values = distributed_utils.GatherFieldToRank0(
                    self.model_part, variable_name, data_location)
                if is_writing_rank:
                    arrays[f"{variable_name}__{data_location}"] = values
            if not is_writing_rank:
                return
        else:
            for variable_name, data_location in self.field_specs:
                variable = Kratos.KratosGlobals.GetVariable(variable_name)
                tensor_adaptor = GetTensorAdaptor(self.model_part, data_location, variable)
                arrays[f"{variable_name}__{data_location}"] = numpy.array(tensor_adaptor.data)

        numpy.savez(self.output_path / f"{self.file_prefix}_{step}.npz", **arrays)
