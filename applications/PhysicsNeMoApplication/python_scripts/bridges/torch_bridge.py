"""Bridge between Kratos tensor adaptors and torch tensors.

torch is an optional runtime dependency: it is imported lazily inside the
functions of this module, so importing the module itself (or the application)
never requires torch to be installed.
"""

import KratosMultiphysics as Kratos


def _TryImportTorch():
    try:
        import torch
        return torch
    except ImportError as e:
        raise ImportError(
            "PhysicsNeMoApplication.torch_bridge requires torch, which could not be imported. "
            "Install it with e.g. 'pip install torch'.") from e


def KratosTensorToTorch(tensor_adaptor):
    """Wraps the data of a tensor adaptor as a torch tensor.

    This is zero-copy: the returned tensor aliases the tensor adaptor's
    internal staging buffer (the same buffer exposed as its numpy ``.data``
    view), NOT the Kratos entity storage itself. Call ``CollectData()`` on the
    adaptor before this function to fill the buffer from the model part, and
    ``StoreData()`` after modifying it to write values back.

    Args:
        tensor_adaptor: A collected Kratos tensor adaptor (e.g.
            Kratos.TensorAdaptors.DoubleTensorAdaptor).

    Returns:
        torch.Tensor: A tensor sharing memory with the adaptor's data buffer.
    """
    torch = _TryImportTorch()
    return torch.from_numpy(tensor_adaptor.data)


def TorchToKratosTensor(torch_tensor, tensor_adaptor, store=True):
    """Copies the values of a torch tensor into a tensor adaptor.

    The copy is unavoidable in the general case: the tensor may live on a GPU
    device or be part of an autograd graph, and ``StoreData()`` in any case
    writes element-wise into the (non-contiguous) Kratos entity storage.

    Args:
        torch_tensor: The source torch tensor. Detached and moved to CPU
            automatically if needed.
        tensor_adaptor: The destination Kratos tensor adaptor. Its shape must
            match the tensor's shape.
        store: If True (default), ``StoreData()`` is called after the copy so
            the values reach the underlying Kratos entities.
    """
    _TryImportTorch()
    array = torch_tensor.detach().cpu().numpy()
    if list(array.shape) != list(tensor_adaptor.data.shape):
        raise ValueError(
            f"Shape mismatch: torch tensor has shape {list(array.shape)} but "
            f"the tensor adaptor expects {list(tensor_adaptor.data.shape)}.")
    tensor_adaptor.data[:] = array
    if store:
        tensor_adaptor.StoreData()
