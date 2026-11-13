"""ONNX deployment bridge: exporting models and running ONNX Runtime sessions.

physicsnemo.deploy.onnx handles the torch -> ONNX export
(export_to_onnx_stream); inference runs through an onnxruntime
InferenceSession created once and cached by the caller (upstream's
run_onnx_inference rebuilds the session on every call, which is wasteful
inside a solution loop, so OnnxInferenceProcess uses CreateOrtSession
instead).

torch/physicsnemo/onnxruntime are optional runtime dependencies, imported
lazily inside the helpers only.
"""


def _TryImportOnnxExport():
    try:
        from physicsnemo.deploy.onnx import export_to_onnx_stream
        return export_to_onnx_stream
    except ImportError as e:
        raise ImportError(
            "PhysicsNeMoApplication.onnx_utils requires physicsnemo (and torch) for ONNX "
            "export, which could not be imported. Install it with e.g. "
            "'pip install nvidia-physicsnemo'.") from e


def _TryImportOnnxRuntime():
    try:
        import onnxruntime
        return onnxruntime
    except ImportError as e:
        raise ImportError(
            "PhysicsNeMoApplication.onnx_utils requires onnxruntime, which could not be "
            "imported. Install it with e.g. 'pip install onnxruntime' (CPU) or "
            "'pip install onnxruntime-gpu' (CUDA).") from e


def ParseDevice(device: str):
    """Splits a device string into (is_cuda, device_index).

    Accepts "cpu", "cuda" and "cuda:N". The index defaults to 0 and is
    validated here rather than being silently dropped.
    """
    text = str(device).strip().lower()
    if not text.startswith("cuda"):
        return False, 0
    remainder = text[len("cuda"):]
    if not remainder:
        return True, 0
    if not remainder.startswith(":") or not remainder[1:].isdigit():
        raise ValueError(
            f"Unsupported device \"{device}\". Use \"cpu\", \"cuda\" or \"cuda:N\".")
    return True, int(remainder[1:])


def CreateOrtSession(model, device: str = "cpu", require_device: bool = False):
    """Creates an onnxruntime InferenceSession for a model.

    Args:
        model: ONNX model bytes, or the path of an .onnx file.
        device: "cpu" (default), "cuda" or "cuda:N".
        require_device: If True, raise when CUDA was requested but the
            session did not get it. The default False keeps onnxruntime's
            own behaviour - CPUExecutionProvider stays appended as a
            fallback - but a warning is logged, because the fallback is
            otherwise completely silent: a missing CUDA build, or simply a
            device index that does not exist, both yield a working
            CPU session that looks like success.

    Returns:
        onnxruntime.InferenceSession
    """
    onnxruntime = _TryImportOnnxRuntime()
    is_cuda, device_index = ParseDevice(device)

    if not is_cuda:
        return onnxruntime.InferenceSession(
            str(model) if not isinstance(model, bytes) else model,
            providers=["CPUExecutionProvider"])

    # ORT resolves its CUDA/cuDNN libraries for free only when torch was
    # imported first; OnnxInferenceProcess deliberately creates the session
    # before importing torch (the whole point being that torch is not
    # needed to deploy), so ask ORT to load them itself.
    if hasattr(onnxruntime, "preload_dlls"):
        try:
            onnxruntime.preload_dlls()
        except Exception:  # never let a best-effort preload break inference
            pass

    session = onnxruntime.InferenceSession(
        str(model) if not isinstance(model, bytes) else model,
        providers=["CUDAExecutionProvider", "CPUExecutionProvider"],
        provider_options=[{"device_id": device_index}, {}])

    if "CUDAExecutionProvider" not in session.get_providers():
        message = (f"device \"{device}\" was requested but the session is running on "
                   f"{session.get_providers()}. Either onnxruntime has no CUDA build "
                   "installed ('pip install onnxruntime-gpu'), or device "
                   f"{device_index} does not exist.")
        if require_device:
            raise RuntimeError(f"CreateOrtSession: {message}")
        import KratosMultiphysics as Kratos
        Kratos.Logger.PrintWarning("PhysicsNeMoApplication.onnx_utils", message)

    return session


_ORT_TYPE_TO_NUMPY = {
    "tensor(float)"   : "float32",
    "tensor(double)"  : "float64",
    "tensor(float16)" : "float16",
}


def NumpyDtypeForOrtInput(ort_input) -> str:
    """Maps an ORT session input's type string to a numpy dtype name
    (defaults to float32 for anything unrecognized)."""
    return _ORT_TYPE_TO_NUMPY.get(ort_input.type, "float32")
