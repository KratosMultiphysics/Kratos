"""Pushes each step's fields into a live queue instead of onto disk.

The export half of the streaming path. Everything it pushes into is defined in
`training.streaming_dataset` - this module is only the Kratos Process that
drives it, so that the training-side machinery stays importable without a
solve and this package keeps its rule that every module here has a Factory.

    solve --StreamingDatasetExportProcess--> LiveSampleQueue
                                                   |
                               KratosStreamingDataset --> TrainModel

Same settings and same gather as DatasetExportProcess, minus the paths.
"""

import KratosMultiphysics as Kratos

from KratosMultiphysics.PhysicsNeMoApplication.training.streaming_dataset import (
    GatherSampleArrays, LiveSampleQueue)


def Factory(settings: Kratos.Parameters, model: Kratos.Model):
    if not isinstance(settings, Kratos.Parameters):
        raise TypeError("Expected input shall be a Parameters object, encapsulating a json string")
    return StreamingDatasetExportProcess(model, settings["Parameters"])


class StreamingDatasetExportProcess(Kratos.Process):
    """Pushes each due step's fields into a queue instead of onto disk.

    Same settings and same gather as DatasetExportProcess, minus the paths.
    Attach the queue with SetQueue (or read the one it creates) and hand it
    to CreateStreamingDataset.

    Settings:
        {
            "model_part_name" : "",
            "list_of_fields"  : [ { "variable_name" : "...",
                                    "data_location" : "node_historical" } ],
            "output_interval" : 1,
            "max_queue_size"  : 0,     // 0 = unbounded
            "close_on_finalize" : true // ends the stream at ExecuteFinalize
        }
    """

    def __init__(self, model: Kratos.Model, settings: Kratos.Parameters) -> None:
        super().__init__()

        default_settings = Kratos.Parameters("""{
            "model_part_name"   : "PLEASE_SPECIFY_MODEL_PART_NAME",
            "list_of_fields"    : [
                {
                    "variable_name" : "PLEASE_SPECIFY_VARIABLE_NAME",
                    "data_location" : "node_historical"
                }
            ],
            "output_interval"   : 1,
            "max_queue_size"    : 0,
            "close_on_finalize" : true
        }""")
        settings.ValidateAndAssignDefaults(default_settings)
        for i in range(settings["list_of_fields"].size()):
            settings["list_of_fields"][i].ValidateAndAssignDefaults(
                default_settings["list_of_fields"][0])

        self.model_part = model[settings["model_part_name"].GetString()]
        self.field_specs = [
            (settings["list_of_fields"][i]["variable_name"].GetString(),
             settings["list_of_fields"][i]["data_location"].GetString())
            for i in range(settings["list_of_fields"].size())
        ]
        self.output_interval = settings["output_interval"].GetInt()
        if self.output_interval < 1:
            raise ValueError(
                f"\"output_interval\" must be >= 1 [ output_interval = {self.output_interval} ].")
        self.close_on_finalize = settings["close_on_finalize"].GetBool()
        self.queue = LiveSampleQueue(settings["max_queue_size"].GetInt())

    def SetQueue(self, queue: LiveSampleQueue) -> None:
        """Uses a caller-provided queue instead of the one created here."""
        self.queue = queue

    def ExecuteFinalizeSolutionStep(self) -> None:
        step = self.model_part.ProcessInfo[Kratos.STEP]
        if step % self.output_interval != 0:
            return
        arrays = GatherSampleArrays(self.model_part, self.field_specs)
        if arrays is not None:   # None on non-writing ranks
            self.queue.Push(arrays)

    def ExecuteFinalize(self) -> None:
        if self.close_on_finalize:
            self.queue.Close()
