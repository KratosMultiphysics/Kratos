"""Training samples streamed out of a RUNNING Kratos solve.

Today training data takes a detour through disk: DatasetExportProcess writes
one .npz per step and CreateNpzDataset reads them back. That is the right
default (samples are reusable, shuffleable, inspectable), but it forces the
solve to finish before training starts and it writes files nobody keeps.

This module removes the detour. A process pushes each step's gathered fields
into a queue; an iterable dataset drains it and yields exactly the same
(inputs, targets) items CreateNpzDataset produces, so downstream consumers
cannot tell the two apart - which the tests assert directly by running the
same case both ways.

    solve ──StreamingDatasetExportProcess──▶ LiveSampleQueue
                                                  │
                              KratosStreamingDataset ──▶ TrainModel

Four upstream contracts are encoded here, each verified rather than assumed:

- physicsnemo's `IterableDatasetBase` is an ABC, **not** a
  torch.utils.data.IterableDataset, so a bare subclass is rejected by torch's
  DataLoader ("has no len()"). The dataset below inherits BOTH.
- physicsnemo's own DataLoader unpacks every yielded item as
  `(data, metadata)`, which would silently swallow the targets of an
  `(inputs, targets)` tuple. Setting `yields_batches = True` bypasses
  collation entirely and passes items through verbatim - correct here anyway,
  since a solver step already produces a whole batch.
- `num_workers > 0` duplicates an iterable stream across workers, so the
  dataset refuses it rather than training on each sample twice.
- Output buffers must not be reused across yields (the loader may still be
  reading the previous one), so every emission allocates fresh tensors.

A stream is single-pass: `TrainModel`'s "streaming" block pins one epoch,
because a second pass over an exhausted queue would otherwise look like a
zero-loss epoch.

torch/physicsnemo are optional runtime dependencies - imported lazily.
"""

import collections
import threading

import numpy

import KratosMultiphysics as Kratos

from KratosMultiphysics.PhysicsNeMoApplication.utilities.tensor_adaptor_dataset_utils import GetTensorAdaptor

_QUEUE_CLOSED = object()


def _TryImportTorch():
    try:
        import torch
        return torch
    except ImportError as e:
        raise ImportError(
            "PhysicsNeMoApplication.streaming_dataset requires torch, which could not be "
            "imported. Install it with e.g. 'pip install torch'.") from e


def _TryImportIterableDatasetBase():
    try:
        from physicsnemo.datapipes import IterableDatasetBase
        return IterableDatasetBase
    except ImportError as e:
        raise ImportError(
            "PhysicsNeMoApplication.streaming_dataset requires physicsnemo >= 2.2 "
            "(IterableDatasetBase landed in 2.2), which could not be imported. Install "
            "it with e.g. 'pip install -U nvidia-physicsnemo'.") from e


def GatherSampleArrays(model_part: Kratos.ModelPart, field_specs) -> dict:
    """One step's fields in DatasetExportProcess's on-disk layout.

    Shared by the file exporter and the live stream so the two cannot drift:
    "TIME"/"STEP" as 0-d arrays, then one entry per field keyed
    "<VARIABLE>__<data_location>". Returns None on non-writing ranks of a
    distributed run, mirroring the exporter.
    """
    step = model_part.ProcessInfo[Kratos.STEP]
    arrays = {
        "TIME": numpy.array(model_part.ProcessInfo[Kratos.TIME]),
        "STEP": numpy.array(step),
    }
    if model_part.IsDistributed():
        from KratosMultiphysics.PhysicsNeMoApplication.distributed import distributed_utils
        is_writing_rank = model_part.GetCommunicator().GetDataCommunicator().Rank() == 0
        for variable_name, data_location in field_specs:
            _, values = distributed_utils.GatherFieldToRank0(
                model_part, variable_name, data_location)
            if is_writing_rank:
                arrays[f"{variable_name}__{data_location}"] = values
        if not is_writing_rank:
            return None
    else:
        for variable_name, data_location in field_specs:
            variable = Kratos.KratosGlobals.GetVariable(variable_name)
            tensor_adaptor = GetTensorAdaptor(model_part, data_location, variable)
            arrays[f"{variable_name}__{data_location}"] = numpy.array(tensor_adaptor.data)
    return arrays


class LiveSampleQueue:
    """Bounded hand-off between a running solve and a training loop.

    The solve is the producer and (in-process) the consumer runs after it, so
    this is deliberately simple: a bounded deque plus a close flag. Bounding
    matters because a solve that outruns training would otherwise grow the
    queue without limit; `Push` blocks once `max_size` is reached.

    Closing is explicit and one-way: a consumer cannot otherwise distinguish
    "no sample yet" from "the solve has ended", which is what would make a
    streaming epoch hang.
    """

    def __init__(self, max_size: int = 0):
        if max_size < 0:
            raise ValueError(f"max_size must be >= 0 (0 = unbounded), got {max_size}.")
        self._items = collections.deque()
        self._max_size = max_size
        self._closed = False
        self._condition = threading.Condition()

    def Push(self, sample) -> None:
        """Appends a sample, blocking while the queue is full."""
        with self._condition:
            if self._closed:
                raise RuntimeError("Push after Close: the queue no longer accepts samples.")
            while self._max_size and len(self._items) >= self._max_size:
                self._condition.wait()
            self._items.append(sample)
            self._condition.notify_all()

    def Close(self) -> None:
        """Signals that no further samples will arrive."""
        with self._condition:
            self._closed = True
            self._condition.notify_all()

    @property
    def closed(self) -> bool:
        with self._condition:
            return self._closed

    def __len__(self) -> int:
        with self._condition:
            return len(self._items)

    def Drain(self):
        """Yields samples until the queue is closed AND empty."""
        while True:
            with self._condition:
                while not self._items and not self._closed:
                    self._condition.wait()
                if not self._items:
                    return
                sample = self._items.popleft()
                self._condition.notify_all()
            yield sample


def _AssembleSample(arrays, keys, torch):
    """The (n_entities, total_width) float32 block CreateNpzDataset builds."""
    blocks = []
    for key in keys:
        if key not in arrays:
            raise KeyError(
                f"Key \"{key}\" is not in the streamed sample; available: "
                f"{sorted(arrays)}.")
        block = numpy.asarray(arrays[key], dtype=numpy.float32)
        blocks.append(block.reshape(block.shape[0], -1))
    # a fresh allocation per emission: the loader may still be reading the
    # previous one, and upstream warns that reusing output buffers races
    return torch.from_numpy(numpy.ascontiguousarray(
        numpy.concatenate(blocks, axis=-1)))


def CreateStreamingDataset(queue: LiveSampleQueue, input_keys, output_keys):
    """An iterable dataset over a LiveSampleQueue.

    Yields `(inputs, targets)` float32 tensors in exactly the layout
    CreateNpzDataset produces - same per-entity flattening, same
    concatenation order - so a streamed run and a dumped run are
    interchangeable.

    Args:
        queue: The LiveSampleQueue a StreamingDatasetExportProcess fills.
        input_keys: Field keys ("<VARIABLE>__<location>") forming the input.
        output_keys: Field keys forming the training target.

    Returns:
        A dataset that is both a physicsnemo IterableDatasetBase and a
        torch.utils.data.IterableDataset.
    """
    torch = _TryImportTorch()
    IterableDatasetBase = _TryImportIterableDatasetBase()

    input_keys = list(input_keys)
    output_keys = list(output_keys)
    if not input_keys or not output_keys:
        raise ValueError("Both input_keys and output_keys must be non-empty.")

    class KratosStreamingDataset(IterableDatasetBase, torch.utils.data.IterableDataset):
        """Live solver output as a single-pass training stream."""

        # physicsnemo's loader would otherwise unpack each item as
        # (data, metadata) and silently discard the targets
        yields_batches = True

        def __init__(self):
            self._queue = queue
            self._epoch = 0
            self._generator = None
            self._emitted = 0

        def __iter__(self):
            worker_info = torch.utils.data.get_worker_info()
            if worker_info is not None and worker_info.num_workers > 1:
                raise RuntimeError(
                    "KratosStreamingDataset does not support num_workers > 1: every "
                    "worker would replay the whole stream, training on each sample "
                    "once per worker. Use num_workers=0.")
            for arrays in self._queue.Drain():
                self._emitted += 1
                yield (_AssembleSample(arrays, input_keys, torch),
                       _AssembleSample(arrays, output_keys, torch))

        def set_epoch(self, epoch: int) -> None:
            self._epoch = int(epoch)
            self._emitted = 0

        def set_generator(self, generator) -> None:
            self._generator = generator

        def close(self) -> None:
            self._queue.Close()

        @property
        def emitted(self) -> int:
            """Samples yielded so far - the streaming counterpart of len()."""
            return self._emitted

    return KratosStreamingDataset()
