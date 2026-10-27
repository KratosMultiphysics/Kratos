"""Opt-in NVTX ranges around the deployment hot paths.

NVTX ranges make the gather/forward/scatter phases of every deployment
process visible in NVIDIA Nsight Systems timelines, next to the solver's
own profile. They are disabled by default and cost nothing until
EnableNvtxRanges() is called (model_registry.LoadModel does so when its
"nvtx_ranges" setting is true); even then a range only emits when torch
is importable and a CUDA device is available, so CPU runs stay no-ops.

torch is imported lazily; module import stays ML-free.
"""

_nvtx_enabled = False


def EnableNvtxRanges() -> None:
    """Turns NVTX range emission on process-wide."""
    global _nvtx_enabled
    _nvtx_enabled = True


def DisableNvtxRanges() -> None:
    """Turns NVTX range emission off process-wide."""
    global _nvtx_enabled
    _nvtx_enabled = False


def NvtxRangesEnabled() -> bool:
    return _nvtx_enabled


class NvtxRange:
    """Context manager wrapping a code block in an NVTX range.

    No-op unless ranges are enabled, torch imports and CUDA is available;
    __exit__ only pops when __enter__ actually pushed, so enabling or
    disabling mid-run can never unbalance the range stack.
    """

    def __init__(self, name: str) -> None:
        self.name = name
        self._pushed = False

    def __enter__(self) -> "NvtxRange":
        if _nvtx_enabled:
            try:
                import torch
            except ImportError:
                return self
            if torch.cuda.is_available():
                torch.cuda.nvtx.range_push(self.name)
                self._pushed = True
        return self

    def __exit__(self, exc_type, exc_value, traceback) -> bool:
        if self._pushed:
            import torch
            torch.cuda.nvtx.range_pop()
            self._pushed = False
        return False
