"""Opt-in CuPy acceleration for the array-heavy bridge paths.

numpy is the reference implementation and the default. CuPy is never
selected implicitly: ``"auto"`` resolves to numpy, and a caller has to ask
for the GPU either per call (a ``"backend"`` setting) or process-wide (the
``KRATOS_PHYSICSNEMO_ARRAY_BACKEND`` environment variable, or
``SetDefaultArrayBackend``). That is deliberate - CuPy changes
floating-point reduction order, so a silent switch would change results
under callers who never asked for it.

Where it pays, and where it does not
------------------------------------
CuPy is not a drop-in numpy replacement here. Most of this application's
per-entity cost is Kratos C++ calls and Python-level container iteration,
which a GPU array library does nothing for, and a host-to-device round trip
costs more than the arithmetic on a small mesh. It pays only where the
array work is genuinely large and either quadratic, GEMM-shaped, or able to
stay on the device and be handed to torch without coming back. Sites that
convert therefore pass a ``size_hint``, and fall back to numpy below the
threshold even when CuPy was explicitly requested.

cupy stays an optional dependency: it is imported lazily inside
``_TryImportCuPy``, so importing this module (or the application) never
requires it.
"""

import os

import numpy

_BACKENDS = ("auto", "numpy", "cupy")
_ENVIRONMENT_VARIABLE = "KRATOS_PHYSICSNEMO_ARRAY_BACKEND"

# Below this many array elements a device round trip costs more than the
# arithmetic it saves. Sites may pass their own measured threshold.
DEFAULT_SIZE_THRESHOLD = 100000

_default_backend = None      # None -> consult the environment variable
_availability = None         # tri-state cache: None unprobed, else bool


def _TryImportCuPy():
    """The cupy module, or an actionable ImportError."""
    try:
        import cupy
        return cupy
    except ImportError as e:
        raise ImportError(
            "PhysicsNeMoApplication.array_backend_utils requires cupy, which could not be "
            "imported. Install the build matching your CUDA toolkit with e.g. "
            "'pip install cupy-cuda12x', or use the \"numpy\" backend.") from e


def IsCuPyAvailable() -> bool:
    """True only when cupy imports *and* a CUDA device actually answers.

    Importability alone is not availability: cupy imports cleanly on a
    machine with no driver and only fails at the first allocation, so this
    probes the runtime. The result is cached, making the per-step check on
    a converted path free.
    """
    global _availability
    if _availability is None:
        try:
            import cupy
            _availability = cupy.cuda.runtime.getDeviceCount() > 0
        except Exception:
            # ImportError, CUDARuntimeError, a driver mismatch - all of them
            # mean the same thing here: use numpy.
            _availability = False
    return _availability


def GetDefaultArrayBackend() -> str:
    """The process-wide backend: SetDefaultArrayBackend, else the env var, else "auto"."""
    if _default_backend is not None:
        return _default_backend
    name = os.environ.get(_ENVIRONMENT_VARIABLE, "auto").strip().lower()
    return name if name in _BACKENDS else "auto"


def SetDefaultArrayBackend(backend) -> None:
    """Sets the process-wide backend, or clears the override with None."""
    global _default_backend
    if backend is not None:
        backend = str(backend).strip().lower()
        if backend not in _BACKENDS:
            raise ValueError(
                f"Unsupported array backend \"{backend}\". Use one of {_BACKENDS}.")
    _default_backend = backend


def ResolveArrayModule(backend: str = "auto", size_hint=None, threshold=None):
    """Picks the array module to run a path with.

    Args:
        backend: "auto" (defer to the process-wide default, which is numpy
            unless asked otherwise), "numpy" or "cupy".
        size_hint: Approximate element count the path will work on. When
            given and below the threshold, numpy is used even if CuPy was
            requested - the transfer would cost more than it saves.
        threshold: The crossover for this site, defaulting to
            DEFAULT_SIZE_THRESHOLD. Sites whose crossover was measured
            elsewhere (rom_bridge's dense basis, for one) pass their own.

    Returns:
        (xp, is_cupy): the array module and whether it is CuPy.
    """
    backend = str(backend).strip().lower()
    if backend not in _BACKENDS:
        raise ValueError(f"Unsupported array backend \"{backend}\". Use one of {_BACKENDS}.")
    if backend == "auto":
        backend = GetDefaultArrayBackend()
        if backend == "auto":
            backend = "numpy"
    if backend == "numpy":
        return numpy, False
    if threshold is None:
        threshold = DEFAULT_SIZE_THRESHOLD
    if size_hint is not None and size_hint < threshold:
        return numpy, False
    if not IsCuPyAvailable():
        return numpy, False
    return _TryImportCuPy(), True


def ToHost(array) -> numpy.ndarray:
    """Any array -> numpy, copying off the device when necessary.

    Kratos tensor adaptors write into host memory they own, so every value
    a converted path produces has to come back through here before it can
    reach a model part.
    """
    if type(array).__module__.startswith("cupy"):
        import cupy
        return cupy.asnumpy(array)
    return numpy.asarray(array)


def ToDevice(array, xp):
    """Moves an array to xp's memory space (a no-op when xp is numpy)."""
    if xp is numpy:
        return ToHost(array)
    return xp.asarray(array)


def AsTorchTensor(array, device=None):
    """Hands an array to torch, without a needless device round trip.

    A CuPy array crosses through DLPack, so a value computed on the GPU
    reaches torch as a CUDA tensor with no copy at all - this is what makes
    a converted path worth having, since the alternative is a
    device-to-host copy immediately followed by a host-to-device one. numpy
    arrays take the ordinary zero-copy ``torch.from_numpy`` route.

    Args:
        array: A numpy or CuPy array.
        device: Optional torch device to move the result to. Left alone when
            omitted, which keeps a CuPy input on the GPU it already lives on.
    """
    try:
        import torch
    except ImportError as e:
        raise ImportError(
            "PhysicsNeMoApplication.array_backend_utils.AsTorchTensor requires torch, which "
            "could not be imported. Install it with e.g. 'pip install torch'.") from e

    if type(array).__module__.startswith("cupy"):
        tensor = torch.from_dlpack(array)
    else:
        tensor = torch.from_numpy(numpy.ascontiguousarray(array))
    return tensor if device is None else tensor.to(device)
