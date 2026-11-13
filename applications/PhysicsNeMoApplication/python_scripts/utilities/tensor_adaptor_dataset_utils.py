"""Factory utilities creating tensor adaptors from a data-location string.

Pure Kratos + numpy: this module never imports torch or physicsnemo.
"""

import typing
import KratosMultiphysics as Kratos

SUPPORTED_DATA_LOCATIONS = (
    "node_historical",
    "node_non_historical",
    "element",
    "condition",
    "element_gauss_point",
    "condition_gauss_point",
)


def GetContainerTensorAdaptor(container, data_location: str, variable: typing.Any,
                              process_info=None, collect: bool = True):
    """Constructs the appropriate core tensor adaptor for an explicit
    entity container (e.g. a communicator LocalMesh's ghost-free Nodes).

    Args:
        container: A Kratos Nodes/Elements/Conditions container matching the
            data location.
        data_location: One of SUPPORTED_DATA_LOCATIONS.
        variable: The Kratos variable to read/write.
        process_info: Required for the Gauss-point locations.
        collect: If True (default), CollectData() is called before returning.

    Returns:
        The constructed tensor adaptor.
    """
    if data_location == "node_historical":
        ta = Kratos.TensorAdaptors.HistoricalVariableTensorAdaptor(container, variable)
    elif data_location in ("node_non_historical", "element", "condition"):
        ta = Kratos.TensorAdaptors.VariableTensorAdaptor(container, variable)
    elif data_location in ("element_gauss_point", "condition_gauss_point"):
        if process_info is None:
            raise ValueError(f"Data location \"{data_location}\" requires a process_info.")
        ta = Kratos.TensorAdaptors.GaussPointVariableTensorAdaptor(container, variable, process_info)
    else:
        raise RuntimeError(
            f"Unsupported data location \"{data_location}\". "
            f"Supported locations: {', '.join(SUPPORTED_DATA_LOCATIONS)}.")

    if collect:
        ta.CollectData()
    return ta


def GetTensorAdaptor(model_part: Kratos.ModelPart, data_location: str, variable: typing.Any,
                     collect: bool = True, local_only: bool = False):
    """Constructs the appropriate core tensor adaptor for a data location.

    Args:
        model_part: The model part whose entities the adaptor reads.
        data_location: One of SUPPORTED_DATA_LOCATIONS.
        variable: The Kratos variable to read/write.
        collect: If True (default), CollectData() is called before returning.
        local_only: If True, read the communicator's LocalMesh containers
            (owned entities only) instead of the model part's own containers,
            which on a distributed model part also include ghosts. This is
            what CouplingInterfaceData does, so it is the layout to use when
            the rows must line up with co-simulation interface data.

    Returns:
        The constructed tensor adaptor.
    """
    source = model_part.GetCommunicator().LocalMesh() if local_only else model_part
    if data_location.startswith("node"):
        container = source.Nodes
    elif data_location.startswith("element"):
        container = source.Elements
    elif data_location.startswith("condition"):
        container = source.Conditions
    else:
        raise RuntimeError(
            f"Unsupported data location \"{data_location}\". "
            f"Supported locations: {', '.join(SUPPORTED_DATA_LOCATIONS)}.")
    return GetContainerTensorAdaptor(container, data_location, variable, model_part.ProcessInfo, collect)


def RowsOfIds(container_ids, query_ids):
    """Row indices of ``query_ids`` inside a container's id array.

    The one implementation of the id -> row lookup the bridges share
    (graph scatter, ROM basis permutation, grid sampling): ``argsort`` +
    ``searchsorted`` rather than a ``{id: row}`` dict plus a generator,
    both of which are interpreter-level loops over every entity. The
    container is NOT assumed id-sorted - the argsort makes it correct
    either way. A query id absent from the container raises KeyError
    naming it.

    Args:
        container_ids: (n,) int array of the container's ids in row order
            (e.g. ``numpy.fromiter((node.Id for node in model_part.Nodes), ...)``).
        query_ids: The ids to locate, any shape; the result has the same shape.

    Returns:
        int64 row indices, ``container_ids[result] == query_ids``.
    """
    import numpy

    part_ids = numpy.asarray(container_ids, dtype=numpy.int64).ravel()
    query = numpy.asarray(query_ids, dtype=numpy.int64)
    flat = query.ravel()
    if part_ids.size == 0:
        if flat.size:
            raise KeyError(int(flat[0]))
        return numpy.empty(query.shape, dtype=numpy.int64)
    order = numpy.argsort(part_ids, kind="stable")
    position = numpy.searchsorted(part_ids[order], flat)
    clipped = numpy.minimum(position, part_ids.size - 1)
    missing = part_ids[order][clipped] != flat
    if missing.any():
        raise KeyError(int(flat[missing][0]))
    return order[clipped].reshape(query.shape)
