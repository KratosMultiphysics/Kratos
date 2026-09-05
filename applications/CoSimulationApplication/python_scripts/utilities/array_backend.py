# CoSimulation imports
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools

# Other imports
import numpy as np
import functools

# tri-state cache for the cupy probe: None => not yet probed, False => unusable, module => usable
_cupy_module = None
# human-readable reason for the last failed probe (None if cupy is usable or not yet probed);
# surfaced in the "cupy" fallback warning so a broken CUDA install does not degrade silently
_cupy_failure_reason = None

def _ResolveCupy():
    """Imports cupy and checks that a CUDA device is visible. The result is cached, this
    function never raises (any failure is treated as "cupy is not usable", with the reason
    cached in "_cupy_failure_reason" for diagnostics).
    """
    global _cupy_module, _cupy_failure_reason
    if _cupy_module is None:
        try:
            import cupy
            if cupy.cuda.runtime.getDeviceCount() > 0:
                _cupy_module = cupy
            else:
                _cupy_module = False
                _cupy_failure_reason = "no CUDA device visible (cupy.cuda.runtime.getDeviceCount() == 0)"
        except Exception as err:
            _cupy_module = False
            _cupy_failure_reason = "{}: {}".format(type(err).__name__, err)
    return _cupy_module or None

def GetCupyFailureReason():
    """Returns the cached reason why cupy is not usable (None if usable or not yet probed)."""
    _ResolveCupy()
    return _cupy_failure_reason

def IsCupyAvailable():
    """Returns whether cupy is importable and at least one CUDA device is visible.
    Used by tests to skip the cupy-specific cases when no GPU is present.
    """
    return _ResolveCupy() is not None


class ArrayBackend:
    """Wraps an array module ("xp", i.e. numpy or cupy) together with the host<->device
    transfers and the scipy-compatible shims required by the convergence accelerators.
    """
    def __init__(self, xp, name):
        self.xp = xp
        self.name = name

    def AsArray(self, array):
        """Uploads a host numpy array to this backend (no-op for the numpy backend)."""
        return self.xp.asarray(array)

    def ToNumpy(self, array):
        """Downloads an array of this backend to host numpy (no-op for the numpy backend)."""
        if self.name == "numpy":
            return array
        return self.xp.asnumpy(array)

    def MakeLinearOperator(self, shape, mat_vec):
        if self.name == "numpy":
            import scipy.sparse.linalg as spla
            return spla.LinearOperator(shape, mat_vec)
        else:
            import cupyx.scipy.sparse.linalg as cpla
            return cpla.LinearOperator(shape, mat_vec)

    def Gmres(self, operator, rhs, atol, rtol):
        if self.name == "numpy":
            import scipy.sparse.linalg as spla
            try:
                return spla.gmres(operator, rhs, atol=atol, rtol=rtol)
            except TypeError: # scipy < 1.12 does not know the "rtol" keyword yet
                return spla.gmres(operator, rhs, atol=atol, tol=rtol)
        else:
            import cupyx.scipy.sparse.linalg as cpla
            # current cupyx is keyword-only "rtol" (no "tol"); keep a "tol" fallback for older cupyx
            try:
                return cpla.gmres(operator, rhs, atol=atol, rtol=rtol)
            except TypeError:
                return cpla.gmres(operator, rhs, atol=atol, tol=rtol)

    def SolveTriangular(self, r_matrix, rhs):
        if self.name == "numpy":
            import scipy.linalg as sla
            return sla.solve_triangular(r_matrix, rhs, check_finite=False)
        else:
            import cupyx.scipy.linalg as csla
            # cupyx accepts "check_finite" as well (default False), pass it explicitly for parity
            return csla.solve_triangular(r_matrix, rhs, check_finite=False)


_backend_cache = {}

def _GetOrCreateBackend(name):
    if name not in _backend_cache:
        xp = np if name == "numpy" else _ResolveCupy()
        _backend_cache[name] = ArrayBackend(xp, name)
    return _backend_cache[name]

def GetArrayBackend(name, context="", echo_level=0):
    """Resolves the "backend" setting ("numpy", "cupy" or "auto") to an ArrayBackend.
    Falls back to numpy (with a warning) if "cupy" was explicitly requested but is not usable
    (cupy not installed, or no CUDA device visible); "auto" falls back silently.
    """
    if name == "numpy":
        return _GetOrCreateBackend("numpy")

    if name == "cupy":
        if IsCupyAvailable():
            return _GetOrCreateBackend("cupy")
        reason = GetCupyFailureReason()
        detail = " ({})".format(reason) if reason else " (cupy is not installed, or no CUDA device was found)"
        cs_tools.cs_print_warning(context, 'Backend "cupy" was requested but is not available{}; falling back to "numpy".'.format(detail))
        return _GetOrCreateBackend("numpy")

    if name == "auto":
        if IsCupyAvailable():
            if echo_level > 0:
                cs_tools.cs_print_info(context, 'cupy is available, using it as the array backend.')
            return _GetOrCreateBackend("cupy")
        return _GetOrCreateBackend("numpy")

    raise Exception('Unknown array backend "{}", available options are "numpy", "cupy" and "auto"!'.format(name))


def HostArrayBoundary(method):
    """Decorator for methods (e.g. "UpdateSolution") that receive and return host numpy
    arrays. Uploads ndarray arguments to "self.backend" before the call and downloads the
    returned array back to host numpy afterwards. Zero-cost passthrough when the backend is
    numpy, so it has no effect on the existing (default) behavior.
    """
    @functools.wraps(method)
    def Wrapper(self, *args, **kwargs):
        backend = self.backend
        if backend.name == "numpy":
            return method(self, *args, **kwargs)
        args = [backend.AsArray(arg) if isinstance(arg, np.ndarray) else arg for arg in args]
        kwargs = {key: (backend.AsArray(val) if isinstance(val, np.ndarray) else val) for key, val in kwargs.items()}
        result = method(self, *args, **kwargs)
        return backend.ToNumpy(result)
    return Wrapper
