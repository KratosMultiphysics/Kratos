#import cupy
#import pyamg
import math
import numpy as np
import scipy.sparse as sp
import KratosMultiphysics as KM
import KratosMultiphysics.FluidDynamicsApplication.python_pool as python_pool
from typing import Any, Callable, Optional, Tuple

# Einsum optimization configuration
opt_type = True

# Backend modules to be defined
xp = None
sparse = None
sparse_linalg = None
asnumpy = None

# Configuration flags
USE_CUPY = None
PRECISION = None
_configured = False # Auxiliary flag to ensure that backend is configured once

USE_AMGX = False #TODO: auxiliary solution to activate/deactivate AMGX until we have a final implementation
GraphBasedSA_AMG = None  # Will be set conditionally when GPU backend is selected

def configure(parallel_type : str, precision : str):
    """
    Auxiliary function to be executed once to set the environment
    """

    # Specify the global variables in the module
    global xp, sparse, sparse_linalg, asnumpy
    global USE_CUPY, PRECISION, _configured

    # Check if the environment has been already configured
    if _configured:
        raise RuntimeError("Environment is already configured.")

    # Backend selection
    if parallel_type == "OpenMP" or parallel_type == "open_mp":
        import scipy.sparse.linalg as sp_linalg

        xp = np
        sparse = sp
        sparse_linalg = sp_linalg
        def _asnumpy(x): return np.asarray(x)
        asnumpy = _asnumpy
        USE_CUPY = False

    elif parallel_type == "GPU" or parallel_type == "gpu":
        import cupy as cp
        import cupyx.scipy.sparse as csp
        import cupyx.scipy.sparse.linalg as csp_linalg

        xp = cp
        sparse = csp
        sparse_linalg = csp_linalg
        asnumpy = cp.asnumpy
        USE_CUPY = True

        import KratosMultiphysics.FluidDynamicsApplication.cupy_amg_linear_solver

    elif parallel_type == "MPI" or parallel_type == "mpi":
        raise ValueError("MPI parallelism is not supported yet.")

    else:
        raise ValueError(f"Unknown parallel type '{parallel_type}'.")

    # Set the floating point precision
    if precision == "float64":
        PRECISION = xp.float64
    elif precision == "float32":
        PRECISION = xp.float32
    else:
        raise ValueError(f"Unknown floating point precision value '{precision}'.")

    # Flag to prevent configure to be performed again
    _configured = True

def reset():
    """
    Auxiliary function to reset the environment in very specific situations (e.g., testing)
    """
    
    # Specify the global variables in the module
    global xp, sparse, sparse_linalg, asnumpy
    global USE_CUPY, PRECISION, _configured

    # Backend modules to be defined
    xp = None
    sparse = None
    sparse_linalg = None
    asnumpy = None

    # Configuration flags
    USE_CUPY = None
    PRECISION = None
    _configured = False # Auxiliary flag to ensure that backend is configured once

def _host_scalars(xp, *scalars):
    """
    Read several 0-d device (or host) scalars back to Python floats with ONE blocking
    device->host transfer. Each value is converted to float64 first; float32 -> float64
    is exact, so the values (and their finiteness) are unchanged.
    """
    stacked = xp.stack([xp.asarray(s, dtype=np.float64).reshape(()) for s in scalars])
    return [float(v) for v in asnumpy(stacked)]


def _require_contiguous(*arrays):
    """
    Validate that arrays can be indexed as packed C buffers by a RawKernel.

    CuPy's matmul honours arbitrary strides, but the fused kernels below use flat
    indices, so a strided view would silently produce wrong results.
    """
    for array in arrays:
        flags = getattr(array, "flags", None)
        if flags is not None and not flags.c_contiguous:
            raise ValueError(
                "Fused GPU kernels require C-contiguous arrays, got shape "
                f"{array.shape} with strides {array.strides}.")


_MASS_MATRIX_KERNELS = {}


def _get_mass_matrix_kernel(dtype_name, n_node, n_comp):
    """
    Compile (once per dtype/shape combination) the element-local mass-matrix GPU kernel.

    One thread handles one (element, component) pair and performs the small
    (n_node x n_node) contraction in registers, replacing a cuBLAS batched GEMM
    that would otherwise serialize all elements in a single CUDA block.

    n_node and n_comp are baked into the source so that the division and the
    inner loops are resolved at compile time.
    """
    key = (dtype_name, n_node, n_comp)
    if key not in _MASS_MATRIX_KERNELS:
        import cupy as cp
        source = r"""
extern "C" __global__ void apply_mass_matrix(
    const DTYPE* __restrict__ x,
    DTYPE* __restrict__ y,
    const DTYPE* __restrict__ M,
    int n_elem)
{
    const int idx = blockIdx.x * blockDim.x + threadIdx.x;
    const int e = idx / N_COMP;
    if (e >= n_elem) return;
    const int l = idx - e * N_COMP;
    #pragma unroll
    for (int i = 0; i < N_NODE; ++i) {
        DTYPE value = 0;
        #pragma unroll
        for (int j = 0; j < N_NODE; ++j) {
            value += M[i * N_NODE + j] * x[(e * N_NODE + j) * N_COMP + l];
        }
        y[(e * N_NODE + i) * N_COMP + l] = value;
    }
}
"""
        source = (source.replace("DTYPE", dtype_name)
                        .replace("N_NODE", str(n_node))
                        .replace("N_COMP", str(n_comp)))
        _MASS_MATRIX_KERNELS[key] = cp.RawKernel(source, "apply_mass_matrix")
    return _MASS_MATRIX_KERNELS[key]


def _apply_mass_matrix_gpu(x, mass_matrix, out):
    """Compute out[e,i,l] = sum_j M[i,j] * x[e,j,l] with a fused GPU kernel."""
    _require_contiguous(x, mass_matrix, out)
    n_elem, n_node = x.shape[0], x.shape[1]
    n_comp = x.size // (n_elem * n_node)
    dtype_name = "float" if x.dtype.itemsize == 4 else "double"
    kernel = _get_mass_matrix_kernel(dtype_name, n_node, n_comp)
    n_threads = n_elem * n_comp
    kernel(((n_threads + 255) // 256,), (256,), (x, out, mass_matrix, n_elem))


_MOMENTUM_STABILIZATION_KERNELS = {}


def _get_momentum_stabilization_kernel(dtype_name, n_node, n_dim):
    """
    Compile (once per dtype/shape combination) the fused momentum-stabilization kernel.

    One thread computes one output entry directly:

        out[e,i,l] = sum_{k,j,m} DN[e,i,k] * a[e,j,k] * M[j,m] * (conv[e,m,l] - pi[e,m,l])

    This replaces the subtract + three batched GEMMs and their two scratch buffers
    with a single launch. n_node and n_dim are compile-time constants.
    """
    key = (dtype_name, n_node, n_dim)
    if key not in _MOMENTUM_STABILIZATION_KERNELS:
        import cupy as cp
        source = r"""
extern "C" __global__ void momentum_stabilization(
    const DTYPE* __restrict__ DN,
    const DTYPE* __restrict__ a,
    const DTYPE* __restrict__ conv,
    const DTYPE* __restrict__ pi,
    const DTYPE* __restrict__ M,
    DTYPE* __restrict__ out,
    int n_elem)
{
    const int idx = blockIdx.x * blockDim.x + threadIdx.x;
    const int l = idx % N_DIM;
    const int i = (idx / N_DIM) % N_NODE;
    const int e = idx / (N_DIM * N_NODE);
    if (e >= n_elem) return;
    DTYPE value = 0;
    #pragma unroll
    for (int k = 0; k < N_DIM; ++k) {
        DTYPE s = 0;
        #pragma unroll
        for (int j = 0; j < N_NODE; ++j) {
            DTYPE t = 0;
            #pragma unroll
            for (int m = 0; m < N_NODE; ++m) {
                t += M[j * N_NODE + m] *
                     (conv[(e * N_NODE + m) * N_DIM + l] - pi[(e * N_NODE + m) * N_DIM + l]);
            }
            s += a[(e * N_NODE + j) * N_DIM + k] * t;
        }
        value += DN[(e * N_NODE + i) * N_DIM + k] * s;
    }
    out[(e * N_NODE + i) * N_DIM + l] = value;
}
"""
        source = (source.replace("DTYPE", dtype_name)
                        .replace("N_NODE", str(n_node))
                        .replace("N_DIM", str(n_dim)))
        _MOMENTUM_STABILIZATION_KERNELS[key] = cp.RawKernel(
            source, "momentum_stabilization")
    return _MOMENTUM_STABILIZATION_KERNELS[key]


def _apply_momentum_stabilization_gpu(dn, a_elemental, conv_elemental, pi_elemental,
                                      mass_matrix, out):
    """Fused GPU implementation of the momentum stabilization term."""
    _require_contiguous(dn, a_elemental, conv_elemental, pi_elemental, mass_matrix, out)
    n_elem, n_node, n_dim = conv_elemental.shape
    dtype_name = "float" if conv_elemental.dtype.itemsize == 4 else "double"
    kernel = _get_momentum_stabilization_kernel(dtype_name, n_node, n_dim)
    n_threads = n_elem * n_node * n_dim
    kernel(((n_threads + 255) // 256,), (256,),
           (dn, a_elemental, conv_elemental, pi_elemental, mass_matrix, out, n_elem))


_LAPLACIAN_KERNELS = {}


def _apply_laplacian_gpu(dn, field, out):
    """Fused out[e,i,l] = sum_{j,k} DN[e,i,k] * DN[e,j,k] * field[e,j,l]."""
    _require_contiguous(dn, field, out)
    n_elem, n_node, n_dim = dn.shape
    n_comp = 1 if field.ndim == 2 else field.shape[2]
    dtype_name = "float" if field.dtype.itemsize == 4 else "double"
    key = (dtype_name, n_node, n_dim, n_comp)
    if key not in _LAPLACIAN_KERNELS:
        import cupy as cp
        source = r"""
extern "C" __global__ void laplacian_apply(
    const DTYPE* __restrict__ DN,
    const DTYPE* __restrict__ field,
    DTYPE* __restrict__ out,
    int n_elem)
{
    const int idx = blockIdx.x * blockDim.x + threadIdx.x;
    const int l = idx % N_COMP;
    const int rest = idx / N_COMP;
    const int i = rest % N_NODE;
    const int e = rest / N_NODE;
    if (e >= n_elem) return;
    DTYPE value = 0;
    const DTYPE* dn_i = DN + (e * N_NODE + i) * N_DIM;
    #pragma unroll
    for (int j = 0; j < N_NODE; ++j) {
        const DTYPE* dn_j = DN + (e * N_NODE + j) * N_DIM;
        DTYPE s = 0;
        #pragma unroll
        for (int k = 0; k < N_DIM; ++k) s += dn_i[k] * dn_j[k];
        value += field[(e * N_NODE + j) * N_COMP + l] * s;
    }
    out[(e * N_NODE + i) * N_COMP + l] = value;
}
"""
        source = (source.replace("DTYPE", dtype_name)
                        .replace("N_NODE", str(n_node))
                        .replace("N_DIM", str(n_dim))
                        .replace("N_COMP", str(n_comp)))
        _LAPLACIAN_KERNELS[key] = cp.RawKernel(source, "laplacian_apply")
    n_threads = n_elem * n_node * n_comp
    _LAPLACIAN_KERNELS[key](((n_threads + 255) // 256,), (256,),
                            (dn, field, out, n_elem))


_ELEMENTAL_GRADIENT_KERNELS = {}


def _apply_elemental_gradient_gpu(dn, field, out):
    """Fused grad[e,c,k] = sum_i field[e,i,c] * DN[e,i,k]."""
    _require_contiguous(dn, field, out)
    n_elem, n_node, n_dim = dn.shape
    n_comp = 1 if field.ndim == 2 else field.shape[2]
    dtype_name = "float" if field.dtype.itemsize == 4 else "double"
    key = (dtype_name, n_node, n_dim, n_comp)
    if key not in _ELEMENTAL_GRADIENT_KERNELS:
        import cupy as cp
        source = r"""
extern "C" __global__ void elemental_gradient(
    const DTYPE* __restrict__ DN,
    const DTYPE* __restrict__ field,
    DTYPE* __restrict__ out,
    int n_elem)
{
    const int idx = blockIdx.x * blockDim.x + threadIdx.x;
    const int k = idx % N_DIM;
    const int rest = idx / N_DIM;
    const int c = rest % N_COMP;
    const int e = rest / N_COMP;
    if (e >= n_elem) return;
    DTYPE value = 0;
    #pragma unroll
    for (int i = 0; i < N_NODE; ++i) {
        value += field[(e * N_NODE + i) * N_COMP + c] * DN[(e * N_NODE + i) * N_DIM + k];
    }
    out[(e * N_COMP + c) * N_DIM + k] = value;
}
"""
        source = (source.replace("DTYPE", dtype_name)
                        .replace("N_NODE", str(n_node))
                        .replace("N_DIM", str(n_dim))
                        .replace("N_COMP", str(n_comp)))
        _ELEMENTAL_GRADIENT_KERNELS[key] = cp.RawKernel(source, "elemental_gradient")
    n_threads = n_elem * n_comp * n_dim
    _ELEMENTAL_GRADIENT_KERNELS[key](((n_threads + 255) // 256,), (256,),
                                     (dn, field, out, n_elem))


_CONVECTIVE_OPERATOR_KERNELS = {}


def _apply_convective_operator_gpu(a_elemental, grad_u, out):
    """Fused conv[e,i,j] = sum_m a[e,i,m] * grad[e,j,m]."""
    _require_contiguous(a_elemental, grad_u, out)
    n_elem, n_node, n_dim = a_elemental.shape
    if grad_u.ndim == 3:
        n_comp = grad_u.shape[1]
        grad_stride = n_dim
    else:
        n_comp = 1
        grad_stride = 1
    dtype_name = "float" if a_elemental.dtype.itemsize == 4 else "double"
    key = (dtype_name, n_node, n_dim, n_comp, grad_stride)
    if key not in _CONVECTIVE_OPERATOR_KERNELS:
        import cupy as cp
        source = r"""
extern "C" __global__ void convective_operator(
    const DTYPE* __restrict__ a,
    const DTYPE* __restrict__ grad,
    DTYPE* __restrict__ out,
    int n_elem)
{
    const int idx = blockIdx.x * blockDim.x + threadIdx.x;
    const int j = idx % N_COMP;
    const int rest = idx / N_COMP;
    const int i = rest % N_NODE;
    const int e = rest / N_NODE;
    if (e >= n_elem) return;
    DTYPE value = 0;
    #pragma unroll
    for (int m = 0; m < N_DIM; ++m) {
        value += a[(e * N_NODE + i) * N_DIM + m] * grad[(e * GRAD_STRIDE + j) * N_DIM + m];
    }
    out[(e * N_NODE + i) * N_COMP + j] = value;
}
"""
        source = (source.replace("DTYPE", dtype_name)
                        .replace("N_NODE", str(n_node))
                        .replace("N_DIM", str(n_dim))
                        .replace("N_COMP", str(n_comp))
                        .replace("GRAD_STRIDE", str(grad_stride)))
        _CONVECTIVE_OPERATOR_KERNELS[key] = cp.RawKernel(source, "convective_operator")
    n_threads = n_elem * n_node * n_comp
    _CONVECTIVE_OPERATOR_KERNELS[key](((n_threads + 255) // 256,), (256,),
                                      (a_elemental, grad_u, out, n_elem))


_NODAL_DIVERGENCE_KERNELS = {}


def _apply_elementwise_nodal_divergence_gpu(dn, uel, shape_functions, out):
    """Fused out[e,I] = N[I] * sum_{j,k} DN[e,j,k] * uel[e,j,k]."""
    _require_contiguous(dn, uel, shape_functions, out)
    n_elem, n_node, n_dim = dn.shape
    dtype_name = "float" if dn.dtype.itemsize == 4 else "double"
    key = (dtype_name, n_node, n_dim)
    if key not in _NODAL_DIVERGENCE_KERNELS:
        import cupy as cp
        source = r"""
extern "C" __global__ void nodal_divergence(
    const DTYPE* __restrict__ DN,
    const DTYPE* __restrict__ uel,
    const DTYPE* __restrict__ N,
    DTYPE* __restrict__ out,
    int n_elem)
{
    const int idx = blockIdx.x * blockDim.x + threadIdx.x;
    const int I = idx % N_NODE;
    const int e = idx / N_NODE;
    if (e >= n_elem) return;
    DTYPE div = 0;
    #pragma unroll
    for (int j = 0; j < N_NODE; ++j) {
        #pragma unroll
        for (int k = 0; k < N_DIM; ++k) {
            const int p = (e * N_NODE + j) * N_DIM + k;
            div += DN[p] * uel[p];
        }
    }
    out[e * N_NODE + I] = N[I] * div;
}
"""
        source = (source.replace("DTYPE", dtype_name)
                        .replace("N_NODE", str(n_node))
                        .replace("N_DIM", str(n_dim)))
        _NODAL_DIVERGENCE_KERNELS[key] = cp.RawKernel(source, "nodal_divergence")
    n_threads = n_elem * n_node
    _NODAL_DIVERGENCE_KERNELS[key](((n_threads + 255) // 256,), (256,),
                                   (dn, uel, shape_functions, out, n_elem))


_PRESSURE_PROJECTION_KERNELS = {}


def _apply_pressure_stabilization_projection_gpu(dn, pi_elemental, shape_functions, out):
    """Fused out[e,I] = sum_k DN[e,I,k] * (sum_J N[J] * Pi[e,J,k])."""
    _require_contiguous(dn, pi_elemental, shape_functions, out)
    n_elem, n_node, n_dim = dn.shape
    dtype_name = "float" if dn.dtype.itemsize == 4 else "double"
    key = (dtype_name, n_node, n_dim)
    if key not in _PRESSURE_PROJECTION_KERNELS:
        import cupy as cp
        source = r"""
extern "C" __global__ void pressure_projection(
    const DTYPE* __restrict__ DN,
    const DTYPE* __restrict__ Pi,
    const DTYPE* __restrict__ N,
    DTYPE* __restrict__ out,
    int n_elem)
{
    const int idx = blockIdx.x * blockDim.x + threadIdx.x;
    const int I = idx % N_NODE;
    const int e = idx / N_NODE;
    if (e >= n_elem) return;
    DTYPE value = 0;
    #pragma unroll
    for (int k = 0; k < N_DIM; ++k) {
        DTYPE t = 0;
        #pragma unroll
        for (int J = 0; J < N_NODE; ++J) {
            t += N[J] * Pi[(e * N_NODE + J) * N_DIM + k];
        }
        value += DN[(e * N_NODE + I) * N_DIM + k] * t;
    }
    out[e * N_NODE + I] = value;
}
"""
        source = (source.replace("DTYPE", dtype_name)
                        .replace("N_NODE", str(n_node))
                        .replace("N_DIM", str(n_dim)))
        _PRESSURE_PROJECTION_KERNELS[key] = cp.RawKernel(source, "pressure_projection")
    n_threads = n_elem * n_node
    _PRESSURE_PROJECTION_KERNELS[key](((n_threads + 255) // 256,), (256,),
                                      (dn, pi_elemental, shape_functions, out, n_elem))


_ELEMENTAL_MEAN_KERNELS = {}


def _apply_elemental_mean_gpu(x, out, inv_n_node):
    """Fused out[e,k] = inv_n_node * sum_i x[e,i,k]."""
    _require_contiguous(x, out)
    n_elem, n_node, n_dim = x.shape
    dtype_name = "float" if x.dtype.itemsize == 4 else "double"
    key = (dtype_name, n_node, n_dim)
    if key not in _ELEMENTAL_MEAN_KERNELS:
        import cupy as cp
        source = r"""
extern "C" __global__ void elemental_mean(
    const DTYPE* __restrict__ x,
    DTYPE* __restrict__ out,
    DTYPE inv_n_node,
    int n_elem)
{
    const int idx = blockIdx.x * blockDim.x + threadIdx.x;
    const int k = idx % N_DIM;
    const int e = idx / N_DIM;
    if (e >= n_elem) return;
    DTYPE value = 0;
    #pragma unroll
    for (int i = 0; i < N_NODE; ++i) {
        value += x[(e * N_NODE + i) * N_DIM + k];
    }
    out[e * N_DIM + k] = value * inv_n_node;
}
"""
        source = (source.replace("DTYPE", dtype_name)
                        .replace("N_NODE", str(n_node))
                        .replace("N_DIM", str(n_dim)))
        _ELEMENTAL_MEAN_KERNELS[key] = cp.RawKernel(source, "elemental_mean")
    n_threads = n_elem * n_dim
    _ELEMENTAL_MEAN_KERNELS[key](((n_threads + 255) // 256,), (256,),
                                 (x, out, x.dtype.type(inv_n_node), n_elem))


class CFDUtils:

    def __init__(self):
        if not _configured:
            raise RuntimeError("Backend not set. Call 'configure()' first.")
        self._pool = None

    # Scratch buffer indices within the internal pool (see allocate_temporaries).
    # Only two buffers exist: every kernel needs at most two live temporaries at
    # once (an element-product tensor plus a small per-element/per-node aux), so
    # each kernel just Get()s the shape slice it needs from one of the two.
    TMP_A = 0  # primary product scratch — sized nelem * max(n_in_el * dim, dim * dim)
    TMP_B = 1  # auxiliary scratch       — sized nelem * dim * dim

    def allocate_temporaries(self, nelem, n_in_el, dim):
        """
        Pre-allocate the two shared temporary buffers used by compute methods.

        The buffers live in a BufferPool with exactly two slots; memory cost is
        nelem * (max(n_in_el * dim, dim * dim) + dim * dim) values, versus one
        buffer per scratch shape before. Compute methods acquire their scratch
        through self._pool.Get() (requesting the shape they need, the pool
        returns a front-slice view) and release it via self._pool.Release()
        before returning, so the pool's in-use tracking makes overlapping or
        nested use fail loudly instead of silently corrupting data. Results are always written to the
        caller-provided (or freshly allocated) 'out' arrays and never alias the
        pool storage.

        Parameters
        ----------
        nelem : int
            Number of elements.
        n_in_el : int
            Number of nodes per element.
        dim : int
            Spatial dimension.
        """
        nd = n_in_el * dim
        dim2 = dim * dim

        # Only two scratch slots are ever live simultaneously (TMP_A product
        # tensor + TMP_B auxiliary); every kernel request fits in these two.
        sizes = [
            nelem * max(nd, dim2),  # TMP_A
            nelem * dim2            # TMP_B
        ]
        self._pool = python_pool.BufferPool(sizes, xp, PRECISION)

    def _ensure_temporaries(self, nelem, n_in_el, dim):
        """Lazily allocate temporaries on first use or if dimensions change.

        The pool's own max_sizes are the single source of truth: if they no
        longer match what (nelem, n_in_el, dim) requires, re-allocate.
        """
        dim2 = dim * dim
        required = [nelem * max(n_in_el * dim, dim2), nelem * dim2]
        if self._pool is None or self._pool.max_sizes != required:
            self.allocate_temporaries(nelem, n_in_el, dim)

    def GetShapeFunctionsOnGaussPoints(self, dim: int, integration_order: int):
        if integration_order == 1:
            # Return shape (ngauss, n_in_el) = (1, dim+1), consistent with integration_order == 2
            N = xp.array([xp.ones(dim+1, dtype=PRECISION) / (dim+1)], dtype=PRECISION)
        elif integration_order == 2:
            if dim == 2:
                N = xp.array([
                    [2.0/3.0, 1.0/6.0, 1.0/6.0],
                    [1.0/6.0, 2.0/3.0, 1.0/6.0],
                    [1.0/6.0, 1.0/6.0, 2.0/3.0]
                ], dtype=PRECISION)
            elif dim == 3:
                N = xp.array([
                    [0.5854101966249685, 0.1381966011250105, 0.1381966011250105, 0.1381966011250105],
                    [0.1381966011250105, 0.5854101966249685, 0.1381966011250105, 0.1381966011250105],
                    [0.1381966011250105, 0.1381966011250105, 0.5854101966249685, 0.1381966011250105],
                    [0.1381966011250105, 0.1381966011250105, 0.1381966011250105, 0.5854101966249685]
                ], dtype=PRECISION)
            else:
                raise Exception(f"Dimension {dim} not supported.")
        else:
            raise Exception("Integration order not supported.")

        return N

    def GetElementalMassMatrix(self, dim):
        if dim == 2:
            M_e = xp.array([
                [2.0/12.0, 1.0/12.0, 1.0/12.0],
                [1.0/12.0, 2.0/12.0, 1.0/12.0],
                [1.0/12.0, 1.0/12.0, 2.0/12.0]
            ], dtype=PRECISION)
        elif dim == 3:
            M_e = xp.array([
                [2.0/20.0, 1.0/20.0, 1.0/20.0, 1.0/20.0],
                [1.0/20.0, 2.0/20.0, 1.0/20.0, 1.0/20.0],
                [1.0/20.0, 1.0/20.0, 2.0/20.0, 1.0/20.0],
                [1.0/20.0, 1.0/20.0, 1.0/20.0, 2.0/20.0]
            ], dtype=PRECISION)
        else:
            raise Exception(f"Dimension {dim} not supported.")

        return M_e

    def GetWallMassMatrix(self, dim):
        if dim == 2:
            M_e = xp.array([
                [1.0/3.0, 1.0/6.0],
                [1.0/6.0, 1.0/3.0]
            ], dtype=PRECISION)
        elif dim == 3:
            M_e = xp.array([
                [2.0/12.0, 1.0/12.0, 1.0/12.0],
                [1.0/12.0, 2.0/12.0, 1.0/12.0],
                [1.0/12.0, 1.0/12.0, 2.0/12.0]
            ], dtype=PRECISION)
        else:
            raise Exception(f"Dimension {dim} not supported.")

        return M_e

    def GetGaussIntegrationWeights(self, dim: int, integration_order: int):
        if integration_order == 1:
            if dim == 2:
                w = xp.array([1.0/2.0],dtype=PRECISION)
            elif dim == 3:
                w = xp.array([1.0/6.0],dtype=PRECISION)
            else:
                raise Exception(f"Dimension {dim} not supported.")
        elif integration_order == 2:
            if dim == 2:
                w = xp.array([
                    1.0/6.0,
                    1.0/6.0,
                    1.0/6.0
                ],dtype=PRECISION)
            elif dim == 3:
                w = xp.array([
                    1.0/24.0,
                    1.0/24.0,
                    1.0/24.0,
                    1.0/24.0
                ],dtype=PRECISION)
            else:
                raise Exception("Dimension not supported.")
        else:
            raise Exception("Integration order not supported.")

        return w

    def ComputeElementalDivergence(self, DN: np.ndarray, uel: np.ndarray, out=None):
        """
        Computes the elemental divergence of a vector field u.

        This term is ∇·u.
        Using Einstein notation: out[e] = sum_k sum_l DN[e,k,l]*uel[e,k,l]

        Parameters
        ----------
        DN : ndarray
            Numpy array with shape (Nelem, n_in_el, dim).
        uel : ndarray
            Numpy array with shape (Nelem, n_in_el, dim).
        out : ndarray, optional
            Pre-allocated output array of shape (Nelem,).
            If None, a new array is allocated.
        """
        # Original einsum (no `out` support in CuPy):
        # return xp.einsum("nij,nij->n", DN, uel, optimize=opt_type)
        nelem, n_in_el, dim = DN.shape
        self._ensure_temporaries(nelem, n_in_el, dim)
        nd = n_in_el * dim

        if out is None:
            out = xp.empty(nelem, dtype=DN.dtype)

        # Zero-copy reshape to 2D
        DN_flat = DN.reshape(nelem, nd)
        uel_flat = uel.reshape(nelem, nd)

        # Element-wise product into scratch buffer, sum rows into caller-provided output
        t_end = self._pool.Get(self.TMP_A, (nelem, nd))
        xp.multiply(DN_flat, uel_flat, out=t_end)
        xp.sum(t_end, axis=1, out=out)
        self._pool.Release(self.TMP_A)

        return out

    def ComputeElementwiseNodalDivergence(self, N: np.ndarray, DN: np.ndarray, uel: np.ndarray, out=None):
        """
        Computes the nodal weighted divergence of a vector field u.

        This term is (w, ∇·u).
        Using Einstein notation: out[e,i] = sum_k sum_l N[i]*DN[e,k,l]*uel[e,k,l]
                                          = xp.einsum("I,eJk,eJk->eI",N,DN,uel,optimize=opt_type)
        Parameters
        ----------
        N : ndarray
            Numpy array with shape (n_in_el,).
        DN : ndarray
            Numpy array with shape (Nelem, n_in_el, dim).
        uel : ndarray
            Numpy array with shape (Nelem, n_in_el, dim).
        out : ndarray, optional
            Pre-allocated output array of shape (Nelem, n_in_el).
            If None, a new array is allocated.
        """
        # Original einsum:
        # return xp.einsum("I,eJk,eJk->eI", N, DN, uel, optimize=opt_type)
        n_in_el = DN.shape[1]
        if N.shape[0] != n_in_el or N.ndim != 1:
            raise Exception("wrong size of N")
        if uel.shape != DN.shape:
            raise Exception("wrong size of uel")

        nelem, n_in_el, dim = DN.shape

        if out is None:
            out = xp.empty((nelem, n_in_el), dtype=DN.dtype)

        # Original implementation (kept for reference / CPU fallback):
        # self._ensure_temporaries(nelem, n_in_el, dim)
        # nd = n_in_el * dim
        # DN_flat = DN.reshape(nelem, nd)
        # uel_flat = uel.reshape(nelem, nd)
        # t_end = self._pool.Get(self.TMP_A, (nelem, nd))
        # t_e = self._pool.Get(self.TMP_B, (nelem,))
        # xp.multiply(DN_flat, uel_flat, out=t_end)
        # xp.sum(t_end, axis=1, out=t_e)
        # xp.multiply(t_e[:, None], N[None, :], out=out)
        # self._pool.Release(self.TMP_A)
        # self._pool.Release(self.TMP_B)

        if USE_CUPY:
            _apply_elementwise_nodal_divergence_gpu(DN, uel, N, out)
        else:
            self._ensure_temporaries(nelem, n_in_el, dim)
            nd = n_in_el * dim
            DN_flat = DN.reshape(nelem, nd)
            uel_flat = uel.reshape(nelem, nd)
            t_end = self._pool.Get(self.TMP_A, (nelem, nd))
            t_e = self._pool.Get(self.TMP_B, (nelem,))
            xp.multiply(DN_flat, uel_flat, out=t_end)
            xp.sum(t_end, axis=1, out=t_e)
            xp.multiply(t_e[:, None], N[None, :], out=out)
            self._pool.Release(self.TMP_A)
            self._pool.Release(self.TMP_B)

        return out

    def ComputeElementalConvectiveOperator(self, a_elemental, grad_u, out=None):
        """
        Computes the convective operator (a·∇u).

        conv[e, i, j] = sum_m a[e,i,m] * DN[e,j,m]
                      = xp.einsum('ekl,eil->eik', grad_u, a_elemental, optimize=opt_type) #vector case
                      = xp.einsum('el,eil->ei', grad_u, a_elemental, optimize=opt_type) #scalar case
        """
        nelem = a_elemental.shape[0]
        n_in_el = a_elemental.shape[1]

        if grad_u.ndim == 3: # vector field
            dim = grad_u.shape[1]
            if out is None:
                out = xp.empty((nelem, n_in_el, dim), dtype=PRECISION)

            # (E, I, L) @ (E, L, K) -> (E, I, K)
            # Original implementation (kept for reference / CPU fallback):
            # xp.matmul(a_elemental, grad_u.swapaxes(1, 2), out=out)
            if USE_CUPY:
                _apply_convective_operator_gpu(a_elemental, grad_u, out)
            else:
                xp.matmul(a_elemental, grad_u.swapaxes(1, 2), out=out)
            return out

        elif grad_u.ndim == 2: # scalar field
            if out is None:
                out = xp.empty((nelem, n_in_el), dtype=PRECISION)

            # (E, I, L) @ (E, L, 1) -> (E, I, 1) -> squeeze to (E, I)
            # Original implementation (kept for reference / CPU fallback):
            # xp.matmul(a_elemental, grad_u[:, :, None], out=out[:, :, None])
            if USE_CUPY:
                _apply_convective_operator_gpu(a_elemental, grad_u, out)
            else:
                xp.matmul(a_elemental, grad_u[:, :, None], out=out[:, :, None])
            return out

        else:
            raise ValueError("grad_u must have 2 dims (scalar) or 3 dims (vector)")

    def ComputeElementalMean(self, v_elemental, out=None):
        """Compute out[e,k] = (1/nnode) * sum_i v_elemental[e,i,k]."""
        nelem, nnode, ndim = v_elemental.shape
        if out is None:
            out = xp.empty((nelem, ndim), dtype=v_elemental.dtype)

        if USE_CUPY:
            _apply_elemental_mean_gpu(v_elemental, out, 1.0 / nnode)
        else:
            # Original implementation (kept for reference / CPU fallback):
            # xp.sum(v_elemental, axis=1, out=out)
            # out /= nnode
            xp.sum(v_elemental, axis=1, out=out)
            out /= nnode

        return out

    def ComputeNDN(self, N: np.ndarray, DN: np.ndarray, pel: np.ndarray, out=None):
        """
        Computes the term (w, ∇p).

        Using Einstein notation: out[e,I,k] = N[I]*DN[e,J,k]*pel[e,J]
                                            = xp.einsum("I,eJk,eJ->eIk", N, DN, pel, optimize=opt_type)

        Parameters
        ----------
        N : ndarray
            Numpy array with shape (n_in_el,).
        DN : ndarray
            Numpy array with shape (Nelem, n_in_el, dim).
        pel : ndarray
            Numpy array with shape (Nelem, n_in_el), i.e., one scalar per node of the element.
        out : ndarray, optional
            Pre-allocated output of shape (Nelem, n_in_el, dim).
        """
        nelem = DN.shape[0]
        n_in_el = DN.shape[1]
        dim = DN.shape[2]

        if pel.shape != (nelem, n_in_el):
            raise ValueError("pel must have shape (nelem, n_in_el) for scalar case. Current shape is:",pel.shape)

        self._ensure_temporaries(nelem, n_in_el, dim)

        if out is None:
            out = xp.empty((nelem, n_in_el, dim), dtype=DN.dtype)

        # Step 1: grad_p[e,k] = sum_J DN[e,J,k] * pel[e,J]
        t_ed = self._pool.Get(self.TMP_B, (nelem, dim))
        xp.matmul(DN.swapaxes(1, 2), pel[:, :, None], out=t_ed[:, :, None])

        # Step 2: out[e,I,k] = N[I] * grad_p[e,k]
        xp.multiply(N[None, :, None], t_ed[:, None, :], out=out)
        self._pool.Release(self.TMP_B)

        return out

    def ComputeDNN(self, N: np.ndarray, DN: np.ndarray, pel: np.ndarray, out=None):
        """
        Computes the term (∇·w, p).

        Using Einstein notation: out[e,I,k] = sum_J DN[e,I,k]*pel[e,J]*N[J]
                                            = xp.einsum("eIk,J,eJ->eIk", DN, N, pel, optimize=opt_type)
        Parameters
        ----------
        N : ndarray
            Numpy array with shape (n_in_el,).
        DN : ndarray
            Numpy array with shape (Nelem, n_in_el, dim).
        pel : ndarray
            Numpy array with shape (Nelem, n_in_el), i.e., one scalar per node of the element.
        out : ndarray, optional
            Pre-allocated output of shape (Nelem, n_in_el, dim).
        """
        # Original einsum:
        # return xp.einsum("eIk,J,eJ->eIk", DN, N, pel, optimize=opt_type)
        nelem, n_in_el, dim = DN.shape
        self._ensure_temporaries(nelem, n_in_el, dim)

        if out is None:
            out = xp.empty((nelem, n_in_el, dim), dtype=DN.dtype)

        # Step 1: s[e] = sum_J pel[e,J] * N[J]  — dot product per element
        t_e = self._pool.Get(self.TMP_B, (nelem,))
        xp.dot(pel, N, out=t_e)

        # Step 2: out[e,I,k] = DN[e,I,k] * s[e]
        xp.multiply(DN, t_e[:, None, None], out=out)
        self._pool.Release(self.TMP_B)

        return out

    def ComputeLaplacianMatrix(self, DN: np.ndarray, out=None):
        """
        Computes the Laplacian local matrix term (∇N, ∇N).

        Using Einstein notation: out[e,i,j] = sum_m DN[e,i,m]*DN[e,j,m]
                                            = xp.einsum("eim,ejm->eij", DN, DN, optimize=opt_type)
        Parameters
        ----------
        DN : ndarray
            Numpy array with shape (Nelem, n_in_el, dim).
        out : ndarray, optional
            Pre-allocated output array, expected to have shape (Nelem, n_in_el, n_in_el).

        Returns
        -------
        out : ndarray
            Laplacian local matrices of shape (Nelem, n_in_el, n_in_el).
        """
        if out is None:
            nelem = DN.shape[0]
            n_in_el = DN.shape[1]
            out = xp.empty((nelem, n_in_el, n_in_el), dtype=DN.dtype)

        # DN (E, n_in_el, dim) @ DN^T (E, dim, n_in_el) -> (E, n_in_el, n_in_el)
        xp.matmul(DN, DN.swapaxes(1, 2), out=out)
        return out

    def ApplyLaplacian(self, DN: np.ndarray, field: np.ndarray, out=None):
        """
        Compute the Laplacian term in a matrix-free manner:
            (∇q, ∇field)   or   (∇w, ∇field)

        Instead of explicitly forming the elemental Laplacian matrix
            L_ij = (∇N_i · ∇N_j),
        which would involve O(nnode²) operations and large temporaries,
        this implementation exploits the factorization:

            L · u = (DN^T · DN) · u
                = DN^T · (DN · u)

        i.e.:
            1) Compute the gradient of the field:
                grad = DN · field
            2) Apply the divergence (projection back to nodes):
                out = DN^T · grad

        This reduces computational cost to O(nnode · dim) and avoids
        constructing the dense (nnode × nnode) Laplacian matrix.

        Parameters
        ----------
        DN : (E, nnode, dim)
            Shape function gradients.
        field : (E, nnode) or (E, nnode, dim)
            Scalar or vector field at element nodes.
        out : ndarray, optional
            Pre-allocated output array.

        Returns
        -------
        (E, nnode) or (E, nnode, dim)
            Laplacian applied to the field.
        """

        nelem, nnode, dim = DN.shape

        if field.ndim == 2:
            out_shape = (nelem, nnode)
        elif field.ndim == 3:
            out_shape = (nelem, nnode, field.shape[2])
        else:
            raise ValueError("Field must have shape (E,nnode) or (E,nnode,dim)")

        if out is None:
            out = xp.empty(out_shape, dtype=field.dtype)

        # Original implementation (kept for reference / CPU fallback):
        # self._ensure_temporaries(nelem, nnode, dim)
        # t_eik = self._pool.Get(self.TMP_A, (nelem, nnode, dim))
        # t_ed = self._pool.Get(self.TMP_B, (nelem, dim))
        # xp.multiply(DN, field[:, :, None], out=t_eik)
        # xp.sum(t_eik, axis=1, out=t_ed)
        # xp.matmul(DN, t_ed[:, :, None], out=out[:, :, None])
        # self._pool.Release(self.TMP_A)
        # self._pool.Release(self.TMP_B)
        # t_edd = self._pool.Get(self.TMP_A, (nelem, dim, dim))
        # grad_view = t_edd[:, :, :ncomp]
        # xp.matmul(DN.swapaxes(1, 2), field, out=grad_view)
        # xp.matmul(DN, grad_view, out=out)
        # self._pool.Release(self.TMP_A)

        if USE_CUPY:
            _apply_laplacian_gpu(DN, field, out)
            return out

        self._ensure_temporaries(nelem, nnode, dim)
        if field.ndim == 2:
            t_eik = self._pool.Get(self.TMP_A, (nelem, nnode, dim))
            t_ed = self._pool.Get(self.TMP_B, (nelem, dim))
            xp.multiply(DN, field[:, :, None], out=t_eik)
            xp.sum(t_eik, axis=1, out=t_ed)
            xp.matmul(DN, t_ed[:, :, None], out=out[:, :, None])
            self._pool.Release(self.TMP_A)
            self._pool.Release(self.TMP_B)
            return out

        ncomp = field.shape[2]
        t_edd = self._pool.Get(self.TMP_A, (nelem, dim, dim))
        grad_view = t_edd[:, :, :ncomp]
        xp.matmul(DN.swapaxes(1, 2), field, out=grad_view)
        xp.matmul(DN, grad_view, out=out)
        self._pool.Release(self.TMP_A)
        return out

    def InterpolateValue(self, N, field, out=None):
        """
        Interpolates field values at Gauss points.
           out =  xp.einsum("gn,en->eg", N, field, optimize=opt_type)

        Parameters
        ----------
        N : (nnode,) or (ngauss, nnode)
        field : (nelem, nnode) or (nelem, nnode, dim)
        out : ndarray, optional
            Pre-allocated output array of shape (nelem, ngauss) or (nelem, ngauss, dim).

        Returns
        -------
        scalar:
            (nelem, ngauss)
        vector:
            (nelem, ngauss, dim)
        """

        # Ensure Gauss dimension exists
        if N.ndim == 1:
            N = N[None, :]   # (1, nnode)

        ngauss = N.shape[0]
        nelem = field.shape[0]

        # Scalar field
        if field.ndim == 2:
            if out is None:
                out = xp.empty((nelem, ngauss), dtype=field.dtype)
            # (nelem, nnode) @ (nnode, ngauss) -> (nelem, ngauss)
            xp.matmul(field, N.T, out=out)
            return out

        # Vector field
        elif field.ndim == 3:
            dim = field.shape[2]
            if out is None:
                out = xp.empty((nelem, ngauss, dim), dtype=field.dtype)
            # N (ngauss, nnode) @ field (nelem, nnode, dim) -> (nelem, ngauss, dim)
            xp.matmul(N, field, out=out)
            return out

        else:
            raise ValueError("Field must have shape (nelem,nnode) or (nelem,nnode,dim)")

    def AssembleVector(self, conn: np.ndarray, vals: np.ndarray, out: np.ndarray):
        """
        Assembles a vector (or block vector) from connectivity and values.

        Parameters
        ----------
        conn : (nelem, nnode)
            Element connectivity.
        vals : (nelem, nnode) or (nelem, nnode, ncomp)
            Element values to assemble.
        out : preallocated array
            Shape must be (conn.max()+1,) or (conn.max()+1, ncomp).
        """
        xp.add.at(out, conn, vals)

    def ComputeElementalGradient(self, DN: np.ndarray, field: np.ndarray, out=None):
        """
        Compute the gradient of a scalar or vector field using element-dependent DN.
        that is for every element we compute:
            grad[k]   = sum_I DN_I/Dx_k p_I  - scalar case = xp.einsum('ei,eik->ek', field, DN, optimize=opt_type)
            grad[k,l] = sum_I DN_I/Dx_l v_Ik - vector case = xp.einsum('eik,eil->ekl', field, DN, optimize=opt_type)
        Parameters
        ----------
        DN : (nelem, nnode, ndim)
            Shape function derivatives for each element.

        field : array
            Field values at element nodes. Must have shape:
            - (nelem, nnode)          for a scalar field
            - (nelem, nnode, ncomp)   for a vector field

        out : ndarray, optional
            Pre-allocated output array of shape (nelem, ndim) or (nelem, ncomp, ndim).

        Notes
        -----
        Computes:
            Scalar: sum_i field[e,i]      * DN[e,i,k]
            Vector: sum_i field[e,i,j]    * DN[e,i,k]

        """

        # ------------------------------
        # Validate DN
        # ------------------------------
        if DN.ndim != 3:
            raise ValueError("DN must have shape (nelem, nnode, ndim)")

        nelem, nnode, ndim = DN.shape

        # ------------------------------
        # Scalar field
        # ------------------------------
        if field.ndim == 2:
            if field.shape != (nelem, nnode):
                raise ValueError("field must have shape (nelem, nnode) for scalar case. Current shape is:",field.shape)

            if out is None:
                out = xp.empty((nelem, ndim), dtype=field.dtype)

            # (E, 1, nnode) @ (E, nnode, ndim) -> (E, 1, ndim) -> squeeze to (E, ndim)
            # Original implementation (kept for reference / CPU fallback):
            # xp.matmul(field[:, None, :], DN, out=out[:, None, :])
            if USE_CUPY:
                _apply_elemental_gradient_gpu(DN, field, out)
            else:
                xp.matmul(field[:, None, :], DN, out=out[:, None, :])
            return out

        # ------------------------------
        # Vector field
        # ------------------------------
        elif field.ndim == 3:
            ncomp = field.shape[2]
            if field.shape[0] != nelem or field.shape[1] != nnode:
                raise ValueError("field must have shape (nelem, nnode, ncomp). Current shape is:",field.shape)

            if out is None:
                out = xp.empty((nelem, ncomp, ndim), dtype=field.dtype)

            # field^T @ DN: (E, ncomp, nnode) @ (E, nnode, ndim) -> (E, ncomp, ndim)
            # Original implementation (kept for reference / CPU fallback):
            # xp.matmul(field.swapaxes(1, 2), DN, out=out)
            if USE_CUPY:
                _apply_elemental_gradient_gpu(DN, field, out)
            else:
                xp.matmul(field.swapaxes(1, 2), DN, out=out)
            return out

        # ------------------------------
        # Invalid field
        # ------------------------------
        raise ValueError("field must have 2 dims (scalar) or 3 dims (vector), Current shape of field is:",field.shape)

    def ComputeBodyForceContribution(self, b_elemental, out=None):
        """
        Compute (w, b) contribution for multiple Gauss points.
                       = xp.einsum("ij,ejd->eid", M_e, b_elemental, optimize=opt_type)
        Parameters
        ----------
        b_elemental : (nelem, nnode, dim)
            body force at element nodes
        out : ndarray, optional
            Pre-allocated output array of shape (nelem, nnode, dim).

        Returns
        -------
        (nelem, ngauss, dim)
        """

        if out is None:
            out = xp.empty((b_elemental.shape), dtype=b_elemental.dtype)

        M_e = self.GetElementalMassMatrix(b_elemental.shape[2])  # (node, node)

        # Original implementation (kept for reference / CPU fallback):
        # xp.matmul(M_e, b_elemental, out=out)   # (e, node, dim)
        if USE_CUPY:
            _apply_mass_matrix_gpu(b_elemental, M_e, out)
        else:
            xp.matmul(M_e, b_elemental, out=out)   # CPU fallback

        return out

    def ComputeConvectiveContribution(self, elem_conv, out=None):
        """
        Compute the elemental convective contribution using a factorized formulation
        based on a nodal convective operator.

        This implementation evaluates the term

            (w, a · ∇u) = ∫_K N_i (a · ∇u) dK

        by exploiting the following assumptions:

        Assumptions
        -----------
        1. Linear finite elements (P1):
        - Shape function gradients ∇N_j are constant within each element.
        - Therefore, the field gradient ∇u is also constant per element:
                ∇u = Σ_j u_j ∇N_j

        2. Finite element interpolation of the convective velocity:
                a(x) = Σ_j N_j(x) a_j

        3. The convective operator is evaluated at nodes:
                c_j = a_j · ∇u

        which implies the FE reconstruction:
                a(x) · ∇u = Σ_j N_j(x) c_j

        Formulation
        -----------
        Using the above, the weak form becomes:

            ∫_K N_i (a · ∇u)
            = ∫_K N_i Σ_j N_j c_j
            = Σ_j (∫_K N_i N_j) c_j
            = Σ_j M_ij c_j

        where M_ij is the consistent elemental mass matrix.

        Therefore, the convective contribution is computed as:

            R_i = Σ_j M_ij c_j

        which is equivalent to standard Galerkin integration under the stated assumptions.

        Parameters
        ----------
        elem_conv : array
            Nodal convective operator per element:

            - Scalar case: shape (n_elem, n_node)
                elem_conv[e, j] = a_j · ∇u

            - Vector case: shape (n_elem, n_node, n_comp)
                elem_conv[e, j, k] = a_j · ∇u_k

        out : ndarray, optional
            Pre-allocated output array. If None, a new array is allocated.
            - Scalar case: shape (n_elem, n_node)
            - Vector case: shape (n_elem, n_node, n_comp)

        Returns
        -------
        array
            Elemental residual contribution:

            - Scalar case: shape (n_elem, n_node)
            - Vector case: shape (n_elem, n_node, n_comp)

            Note:
            The returned values correspond to the reference element contribution.
            Any geometric scaling (e.g. detJ or element volume) must be applied externally.

        Notes
        -----
        - This formulation avoids explicit quadrature over Gauss points and replaces it
        with a matrix-vector (or matrix-matrix) product using the elemental mass matrix.

        - It is algebraically equivalent to standard Gauss integration if:
            * the same FE interpolation is used for the velocity field,
            * the gradient is computed consistently,
            * and the mass matrix corresponds to the same reference element.

        - Differences with Gauss-based implementations can arise if:
            * under-integration is used,
            * the velocity field is not interpolated from nodal values,
            * or inconsistent geometric scaling is applied.

        - This formulation is particularly well-suited for GPU execution since it:
            * avoids per-Gauss-point loops,
            * uses small dense tensor contractions,
            * and enables further fusion with gradient computation.
        """

        if elem_conv.ndim == 3: # vector field
            nelem, n_node, n_dim = elem_conv.shape
            M_e = self.GetElementalMassMatrix(n_dim) #elemental mass matrix (n_node, n_node)

            if out is None:
                out = xp.empty((nelem, n_node, n_dim), dtype=elem_conv.dtype)
            elif out.shape != (nelem, n_node, n_dim):
                raise ValueError("out must have shape (n_elem, n_node, n_comp)")

            # Original implementation (kept for reference / CPU fallback):
            # xp.matmul(M_e, elem_conv, out=out)
            if USE_CUPY:
                _apply_mass_matrix_gpu(elem_conv, M_e, out)
            else:
                xp.matmul(M_e, elem_conv, out=out)   # CPU fallback
            return out

        elif elem_conv.ndim == 2: # scalar field
            nelem, n_node = elem_conv.shape
            dim = n_node - 1
            M_e = self.GetElementalMassMatrix(dim) #elemental mass matrix (n_node, n_node)

            if out is None:
                out = xp.empty((nelem, n_node), dtype=elem_conv.dtype)
            elif out.shape != (nelem, n_node):
                raise ValueError("out must have shape (n_elem, n_node)")

            # Original implementation (kept for reference / CPU fallback):
            # xp.matmul(elem_conv, M_e.T, out=out)
            if USE_CUPY:
                # The kernel computes M @ x; since M_e.T[j,i] == M_e[i,j],
                # this is equivalent to the original elem_conv @ M_e.T.
                _apply_mass_matrix_gpu(elem_conv, M_e, out)
            else:
                xp.matmul(elem_conv, M_e.T, out=out)   # CPU fallback
            return out

        else:
            raise ValueError("Provided elemental nodal convective operator must have 2 dims (scalar) or 3 dims (vector)")

    def ComputeMomentumStabilization(self, DN, a_elemental, conv_elemental, pi_elemental, out=None):
        """
        Compute the momentum stabilization term:

            (a · ∇w, a · ∇u - Π)

        using a consistent mass-matrix projection instead of Gauss quadrature loop.

        The nodal quantity

            β_j = (a_j · ∇u) - Π_j

        is projected with the elemental mass matrix, yielding:

            R_{i,l} = Σ_{j,m,k} DN_{i,k} M_{j m} β_{m,k} a_{j,l}

        Assumptions
        -----------
        - Linear elements (constant ∇u per element)
        - FE interpolation of convective operator a · ∇u and projection Π

        Parameters
        ----------
        DN : (n_elem, n_node, n_dim)
        a_elemental : (n_elem, n_node, n_dim)
        conv_elemental : (n_elem, n_node, n_dim)
        pi_elemental : (n_elem, n_node, n_dim)
        out : ndarray, optional
            Pre-allocated output array of shape (n_elem, n_node, n_dim). If None, a new array is allocated.

        Returns
        -------
        (n_elem, n_node, n_dim)

        Notes
        -----
        - Equivalent to Gauss integration under consistent interpolation.
        - Geometric scaling (detJ / volume) applied externally.
        - Uses batched matrix operations to eliminate loops and temporary 4th-order tensors.
        - Employs out parameter to avoid implicit temporaries in matmul.
        """

        n_elem, n_node, n_dim = conv_elemental.shape

        if out is None:
            out = xp.empty((n_elem, n_node, n_dim), dtype=conv_elemental.dtype)
        elif out.shape != (n_elem, n_node, n_dim):
            raise ValueError(
                f"out must have shape (n_elem, n_node, n_dim) = ({n_elem}, {n_node}, {n_dim}), "
                f"but got shape {out.shape}"
            )
        # Note: out is fully overwritten by the kernel / matmul,
        # so no zeroing is needed when reusing pre-allocated arrays.

        # Original computation (commented for reference):
        # S = xp.zeros((n_elem, n_dim, n_dim), dtype=conv_elemental.dtype)
        # for j in range(n_node):
        #     for m in range(n_node):
        #         w = M_e[j, m]
        #         S += w * (a_elemental[:, j, :, None] * beta[:, m, None, :])
        # out = xp.matmul(DN, S, out=out)

        # Vectorized computation (no loops, no 4th-order tensor):
        # out[e,i,l] = sum_{j,m,k} DN[e,i,k] * M[j,m] * a[e,j,k] * beta[e,m,l]
        # Factor as: S[e,k,l] = sum_{j,m} M[j,m] * a[e,j,k] * beta[e,m,l]
        # Then: out[e,i,l] = sum_k DN[e,i,k] * S[e,k,l]

        # Original implementation (kept for reference / CPU fallback):
        # beta = conv_elemental - pi_elemental
        # n_elem, n_node, n_dim = beta.shape
        # M_e = self.GetElementalMassMatrix(n_dim)
        # self._ensure_temporaries(n_elem, n_node, n_dim)
        # t_eik = self._pool.Get(self.TMP_A, (n_elem, n_node, n_dim))
        # t_edd = self._pool.Get(self.TMP_B, (n_elem, n_dim, n_dim))
        # xp.matmul(M_e[None, :, :], beta, out=t_eik)
        # xp.matmul(a_elemental.swapaxes(1, 2), t_eik, out=t_edd)
        # xp.matmul(DN, t_edd, out=out)
        # self._pool.Release(self.TMP_A)
        # self._pool.Release(self.TMP_B)

        if USE_CUPY:
            M_e = self.GetElementalMassMatrix(n_dim)
            _apply_momentum_stabilization_gpu(
                DN, a_elemental, conv_elemental, pi_elemental, M_e, out)
        else:
            beta = conv_elemental - pi_elemental
            M_e = self.GetElementalMassMatrix(n_dim)
            self._ensure_temporaries(n_elem, n_node, n_dim)
            t_eik = self._pool.Get(self.TMP_A, (n_elem, n_node, n_dim))
            t_edd = self._pool.Get(self.TMP_B, (n_elem, n_dim, n_dim))
            xp.matmul(M_e[None, :, :], beta, out=t_eik)
            xp.matmul(a_elemental.swapaxes(1, 2), t_eik, out=t_edd)
            xp.matmul(DN, t_edd, out=out)
            self._pool.Release(self.TMP_A)
            self._pool.Release(self.TMP_B)

        return out

    def ComputeDivDivStabilization(self, N: np.array, DN: np.ndarray, u_elemental : np.ndarray, Pi_div_elemental: np.ndarray, out=None):
        """
        Compute div-div stabilization term:

            (∇·w, ∇·u - Π)

        Strategy:
            1) Compute divergence of u
            2) Compute projected Π
            3) Form scalar residual
            4) Apply DN once

        Parameters
        ----------
        N : (nnode,)
        DN : (E, nnode, dim)
        u_elemental : (E, nnode, dim)
        Pi_div_elemental : (E, nnode)
        out : ndarray, optional
            Pre-allocated output array of shape (E, nnode, dim).

        Returns
        -------
        (E, nnode, dim)
        """
        nelem, nnode, dim = DN.shape
        self._ensure_temporaries(nelem, nnode, dim)

        if out is None:
            out = xp.empty((nelem, nnode, dim), dtype=DN.dtype)

        # Two-slot layout: TMP_B column 0 holds div(u), column 1 holds div(Π);
        # TMP_A carries the product tensors (released between the two sums).
        t_b = self._pool.Get(self.TMP_B, (nelem, dim * dim))

        # Step 1: divergence of u → (E,)
        t_a = self._pool.Get(self.TMP_A, (nelem, nnode, dim))
        xp.multiply(DN, u_elemental, out=t_a)
        xp.sum(t_a, axis=(1, 2), out=t_b[:, 0])
        self._pool.Release(self.TMP_A)

        # Step 2: projection of Π → (E,)
        t_a = self._pool.Get(self.TMP_A, (nelem, nnode))
        xp.multiply(N[None, :], Pi_div_elemental, out=t_a)
        xp.sum(t_a, axis=1, out=t_b[:, 1])
        self._pool.Release(self.TMP_A)

        # Step 3: scalar residual → (E,) — in-place subtraction
        xp.subtract(t_b[:, 0], t_b[:, 1], out=t_b[:, 0])

        # Step 4: apply DN → (E, nnode, dim)
        xp.multiply(DN, t_b[:, 0][:, None, None], out=out)
        self._pool.Release(self.TMP_B)

        return out

    def ComputePressureStabilizationProjectionTerm(self, N: np.ndarray, DN: np.ndarray, Pi_press_el: np.ndarray, out=None):
        """
        implements (∇q,Pi_pressure)

        using the einsum relation:
            out[e,I] = sum_{J,k} DN[e,I,k] * N[J] * Pi_press_el[e,J,k]
            = xp.einsum("eIk,J,eJk->eI",DN,N,Pi_press_el, optimize=opt_type)

        Parameters
        ----------
        N : (n_node,)
            shape function values
        DN : (nelem, n_node, dim)
            shape function derivatives
        Pi_press_el : (nelem, n_node, dim)
            pressure projection at element nodes
        out : ndarray, optional
            Pre-allocated output array of shape (nelem, n_node).

        Returns
        -------
        (nelem, n_node)
            Scalar result per element-node pair
        """
        nelem, n_node, n_dim = DN.shape

        if out is None:
            out = xp.empty((nelem, n_node), dtype=DN.dtype)

        # Original implementation (kept for reference / CPU fallback):
        # self._ensure_temporaries(nelem, n_node, n_dim)
        # t_ed = self._pool.Get(self.TMP_B, (nelem, n_dim))
        # t_eik = self._pool.Get(self.TMP_A, (nelem, n_node, n_dim))
        # xp.matmul(Pi_press_el.swapaxes(1, 2), N, out=t_ed)
        # xp.multiply(DN, t_ed[:, None, :], out=t_eik)
        # xp.sum(t_eik, axis=2, out=out)
        # self._pool.Release(self.TMP_B)
        # self._pool.Release(self.TMP_A)

        if USE_CUPY:
            _apply_pressure_stabilization_projection_gpu(DN, Pi_press_el, N, out)
        else:
            self._ensure_temporaries(nelem, n_node, n_dim)
            t_ed = self._pool.Get(self.TMP_B, (nelem, n_dim))
            t_eik = self._pool.Get(self.TMP_A, (nelem, n_node, n_dim))
            xp.matmul(Pi_press_el.swapaxes(1, 2), N, out=t_ed)
            xp.multiply(DN, t_ed[:, None, :], out=t_eik)
            xp.sum(t_eik, axis=2, out=out)
            self._pool.Release(self.TMP_B)
            self._pool.Release(self.TMP_A)

        return out

    def ComputeNavierSlipWallContribution(self, v_elemental_wall, elemental_tang_proj, viscosity, slip_length, out=None):
        """
        Compute the contribution of Navier slip wall boundary conditions to the residual.

        The Navier slip condition can be expressed as:

            u_t = L_s * (∂u_t/∂n)
        where u_t is the tangential velocity at the wall, L_s is the slip length, and ∂u_t/∂n is the normal derivative of the tangential velocity.
        The weak form contribution for a test function w is:
            R = ∫_Γ w · (u_t - L_s * (∂u_t/∂n)) dΓ
        This implementation computes the contribution to the residual vector R for nodes on the slip boundary.
        Parameters
        ----------
        v_elemental_wall : (n_node, dim)
            Wall velocity at the boundary nodes.
        elemental_tang_proj : (n_node, dim, dim)
            Tangential projector for the wall elements.
        viscosity : float
            Dynamic viscosity for the Navier slip condition.
        slip_length : float
            Slip length L_s for the Navier slip condition.
        out : ndarray, optional
            Pre-allocated output array of shape (n_node, dim) to store the contribution.
        Returns
        -------
        (n_node, dim)
            Contribution to the residual vector for the slip boundary condition.
        """

        if out is None:
            out = xp.empty((v_elemental_wall.shape), dtype=v_elemental_wall.dtype)

        dim = v_elemental_wall.shape[1]
        M_w = self.GetWallMassMatrix(dim)  # (n_node, n_node)

        v_tang_elem = xp.einsum("eij,enj->eni", elemental_tang_proj, v_elemental_wall, optimize=opt_type)
        out = xp.matmul(M_w, v_tang_elem, out=out)
        out *= -viscosity / (slip_length + 1e-12) # Avoid division by zero in the no-slip limit

        return out

    #REMARK: this function has to be run on the cpu AND NOT ON THE GPU - it will return a gpu matrix if that is the case
    def AllocateScalarMatrix(self,conn : np.ndarray):
        #TODO: make it efficient
        size = xp.max(conn)+1
        Agraph = KM.SparseContiguousRowGraph(size)
        for e in range(conn.shape[0]):
            Agraph.AddEntries(conn[e])
        Agraph.Finalize()

        #assembling matrix
        A = KM.CsrMatrix(Agraph)

        if USE_CUPY:
            assembly_indices = A.GetEquationIdCsrIndices(conn)
            A_cu = sparse.csr_matrix((
                xp.asarray(
                    A.value_data(), dtype=PRECISION),
                xp.asarray(
                    A.index2_data(),dtype=xp.int32),
                xp.asarray(
                    A.index1_data(),dtype=xp.int32)),
                shape=(A.Size1(), A.Size2())
            )
            return A_cu,xp.asarray(assembly_indices,dtype=xp.int32)
        else:
            assembly_indices = A.GetEquationIdCsrIndices(conn)
            return A, xp.asarray(assembly_indices,dtype=xp.int32)

    def AssembleScalarMatrix(self, LHSel, conn, A : KM.CsrMatrix):
        A.SetValue(0.0)
        A.BeginAssemble()
        for e in range(conn.shape[0]):
            A.Assemble(LHSel[e],conn[e])
        A.FinalizeAssemble()

    def AssembleScalarMatrixByCSRIndices(self, LHSel, csr_indices, A):
        if USE_CUPY:
            A.data.fill(0.0)
            self.AssembleVector(csr_indices.ravel(), LHSel.ravel(), A.data)
        else:
            A.value_data().fill(0.0)
            self.AssembleVector(csr_indices.ravel(), LHSel.ravel(), A.value_data())

    def GetScalarMatrixDirichletIndices(self, fixed_values, LHS):
        """
        Returns the indices to apply the sparse matrix BCs
        Note that BCs are applied in a block-based manner
        """

        # Get sparse matrix vectors
        if USE_CUPY:
            n = LHS.shape[0]
            data = LHS.data
            indices = LHS.indices
            indptr = LHS.indptr
        else:
            n = LHS.Size1()
            data = LHS.value_data()
            indices = LHS.index2_data()
            indptr = LHS.index1_data()

        # 1. Create a boolean lookup for fixed indices (O(1) lookup speed)
        is_fixed = xp.zeros(n, dtype=bool)
        is_fixed[fixed_values] = True

        # 2. Map every entry in the sparse 'data' array to its row index
        row_start_flags = xp.zeros(data.size, dtype=bool)
        row_start_flags[indptr[:-1]] = True
        row_map = xp.cumsum(row_start_flags) - 1

        # 3. Create masks for entries in fixed rows and fixed columns
        mask_row = is_fixed[row_map]
        mask_col = is_fixed[indices]

        # 4. Get the indices of the fixed entries (all entries in fixed rows or columns)
        # This is the most expensive step, done in a single vectorized pass
        csr_data_mask = mask_row | mask_col
        csr_data_indices = xp.nonzero(csr_data_mask)[0].astype(xp.int32)

        # 5. Get the indices of the diagonal entries (only the indices in the fixed rows diagonal)
        csr_diag_mask = (row_map == indices) & mask_row
        csr_diag_indices = xp.nonzero(csr_diag_mask)[0].astype(xp.int32)

        return csr_data_indices, csr_diag_indices

    def ApplyHomogeneousDirichlet(self, fixed_values, csr_data_indices, csr_diag_indices, LHS, RHS, diag_value=1.0):
        """
        Applies homogeneous Dirichlet BCs (value=0) to a CSR matrix and RHS.
        Maintains symmetry by zeroing both rows and columns of fixed indices.
        """

        # Get sparse matrix value vector
        if USE_CUPY:
            data = LHS.data
        else:
            data = LHS.value_data()

        # Make zero all fixed value entries
        xp.put(data, csr_data_indices, 0.0)

        # Restore the diagonal for the fixed rows
        xp.put(data, csr_diag_indices, diag_value)

        # Set RHS to 0 at fixed indices
        xp.put(RHS, fixed_values, 0.0)

        return LHS, RHS

    def ConstructPreconditioner(self, A):
        """
        Constructs an algebraic preconditioner for the given matrix A.
        For simplicity, we use the Jacobi (diagonal) preconditioner here.
        """
        if USE_CUPY:
            return KM.FluidDynamicsApplication.cupy_amg_linear_solver.GraphBasedSA_AMG(A, max_coarse=300)
        else:
            return None
    
    def robust_cg(
        self,
        A,
        b,
        x0 = None,
        *,
        xp,
        rtol = 1e-5,
        atol = 0.0,
        maxiter = None,
        M = None,
        callback = None,
        recompute_residual_every: int = 50,
        residual_gap_tol: float = 0.1,
        replace_residual: bool = True,
        restart_on_residual_replacement: bool = True,
        use_mixed_reductions: bool = True,
        return_info_dict: bool = False,
    ):
        """
        Conjugate Gradient with periodic true-residual checks.

        Parameters
        ----------
        A:
            Matrix, sparse matrix, LinearOperator-like object, or callable.
            Must represent a Hermitian/symmetric positive-definite operator.

        b:
            Right-hand side vector.

        x0:
            Initial guess.

        xp:
            Array namespace, typically `cupy` or `numpy`.

        rtol, atol:
            Convergence is accepted when
                ||b - A x|| <= max(atol, rtol * ||b||).
            Note that this is checked using the true residual, not only the
            recursively updated residual.

        maxiter:
            Maximum number of iterations. Defaults to 10 * n.

        M:
            Preconditioner. May be None, callable, matrix-like, or object with
            `.matvec`.

        callback:
            Called as callback(x) after each iteration.

        recompute_residual_every:
            Recompute r_exact = b - A x every this many iterations. Use 0 or None
            to disable periodic checks. Even when disabled, the true residual is
            still computed before accepting convergence from the recursive residual.

        residual_gap_tol:
            If
                ||r_exact - r|| / max(||r_exact||, tiny) > residual_gap_tol,
            the recursive residual is considered unreliable.

        replace_residual:
            Replace r by r_exact when the residual gap is too large.

        restart_on_residual_replacement:
            Restart the CG direction after residual replacement. This is the safer
            choice because replacing r disturbs the strict recurrence.

        use_mixed_reductions:
            If True, use FP32 products/squares with FP64 accumulation for dot/norm.
            If False, use xp.vdot and xp.linalg.norm.

        return_info_dict:
            If False, return SciPy/CuPy-style info integer:
            0 on convergence, iteration count on failure.
            If True, return a dictionary with residual histories and reason.

        Returns
        -------
        x, info
            x is the approximate solution. info is either an integer or a dict.
        """

        def dot_mixed(xp, a, b):
            """
            Dot product for real-valued vectors with FP32 product and FP64 accumulation.

            For complex vectors, replace this with an xp.vdot-based implementation.
            """
            return xp.sum(a * b, dtype=xp.float64)


        def norm_mixed(xp, a):
            """
            Euclidean norm with FP32 square and FP64 accumulation.
            """
            return xp.sqrt(xp.sum(a * a, dtype=xp.float64))


        def dot_default(xp, a, b):
            """
            Default namespace dot product. Useful if you do not want mixed reductions.
            """
            return xp.vdot(a, b)


        def norm_default(xp, a):
            """
            Default namespace norm. Useful if you do not want mixed reductions.
            """
            return xp.linalg.norm(a)


        def make_system(A, M, x0, b, xp):
            """
            Small local replacement for CuPy/SciPy's private _make_system helper.

            It supports:
            - A as an object with .matvec(v)
            - A as a callable A(v)
            - A as an array/sparse matrix supporting A @ v

            Same for M. If M is None, identity preconditioning is used.
            """
            b = xp.asarray(b)

            n = b.shape[0]
            if x0 is None:
                x = xp.zeros_like(b)
            else:
                x = xp.asarray(x0).copy()

            if hasattr(A, "matvec"):
                matvec = A.matvec
                shape = A.shape
            elif callable(A):
                matvec = A
                shape = (n, n)
            else:
                matvec = lambda v: A @ v
                shape = A.shape

            if M is None:
                psolve = lambda r: r
            elif hasattr(M, "matvec"):
                psolve = M.matvec
            elif callable(M):
                psolve = M
            else:
                psolve = lambda r: M @ r

            class Operator:
                def __init__(self, shape, matvec):
                    self.shape = shape
                    self.matvec = matvec

            class Preconditioner:
                def __init__(self, matvec):
                    self.matvec = matvec

            return Operator(shape, matvec), Preconditioner(psolve), x, b





        Aop, Mop, x, b = make_system(A, M, x0, b, xp)
        matvec = Aop.matvec
        psolve = Mop.matvec

        n = Aop.shape[0]
        if maxiter is None:
            maxiter = n * 10

        if n == 0:
            x_empty = xp.empty_like(b)
            if return_info_dict:
                return x_empty, {"converged": True, "iterations": 0, "reason": "empty system"}
            return x_empty, 0

        dot = dot_mixed if use_mixed_reductions else dot_default
        norm = norm_mixed if use_mixed_reductions else norm_default

        b_norm = norm(xp, b)
        b_norm_host = float(b_norm)

        if b_norm_host == 0.0:
            # Returning 'b' itself would alias the caller's buffer (e.g., a pool slot that is
            # released right after the solve and then re-acquired). The exact solution of
            # A x = 0 is x = 0, so return a fresh zero vector instead.
            # Original (kept for reference):
            # if return_info_dict:
            #     return b, {"converged": True, "iterations": 0, "reason": "zero rhs"}
            # return b, 0
            x_zero = xp.zeros_like(b)
            if return_info_dict:
                return x_zero, {"converged": True, "iterations": 0, "reason": "zero rhs"}
            return x_zero, 0

        stop_atol = max(float(atol), float(rtol) * b_norm_host)

        # Initial true residual.
        r = b - matvec(x)
        resid = norm(xp, r)
        resid_host = float(resid)

        recursive_residual_norms = [resid_host]
        true_residual_norms = [resid_host]
        residual_gap_ratios = []
        replacement_iterations = []

        if resid_host <= stop_atol:
            if return_info_dict:
                return x, {
                    "converged": True,
                    "iterations": 0,
                    "reason": "initial guess satisfies tolerance",
                    "residual_norm": resid_host,
                    "recursive_residual_norms": recursive_residual_norms,
                    "true_residual_norms": true_residual_norms,
                    "residual_gap_ratios": residual_gap_ratios,
                    "replacement_iterations": replacement_iterations,
                }
            return x, 0

        iters = 0
        rho = 0
        restart_p = True
        tiny = 1e-300

        reason = "maximum iterations reached"
        converged = False

        while iters < maxiter:
            z = psolve(r)

            rho1 = rho
            rho = dot(xp, r, z)

            # ORIGINAL per-scalar blocking reads and checks (kept for reference):
            # rho_host = float(rho)
            # # if not bool(xp.isfinite(rho)) or rho_host <= 0.0:
            # #     print("norm of r = ",norm(xp,r))
            # #     print("norm of z=M@r= ",norm(xp,z))
            # #     print("bool(xp.isfinite(rho))= ",bool(xp.isfinite(rho)))
            # #     print("rho_host= ",rho_host)
            # #     reason = "non-positive or non-finite r^T M r"
            # #     break
            #
            # if not bool(xp.isfinite(rho)):
            #     print("norm of r = ",norm(xp,r))
            #     print("norm of z=M@r= ",norm(xp,z))
            #     print("bool(xp.isfinite(rho))= ",bool(xp.isfinite(rho)))
            #     print("rho_host= ",rho_host)
            #     reason = "non-finite r^T M r"
            #     break
            #
            #
            # if rho_host < 0.0: 
            #     print("norm of r = ",norm(xp,r))
            #     print("norm of z=M@r= ",norm(xp,z))
            #     print("bool(xp.isfinite(rho))= ",bool(xp.isfinite(rho)))
            #     print("rho_host= ",rho_host)
            #     print("rho_host < 0 - going on as without preconditioner")
            #     reason = " r^T M r < 0"
            #     break
            #
            # #just proceed as nothing in the case it is exactly 0
            # if rho_host == 0.0:
            #     print("rho_host == 0 - simply proceeding even though this should not happen")
            #     pass
            #
            #
            # if restart_p:
            #     p = z.copy() if hasattr(z, "copy") else z
            #     restart_p = False
            # else:
            #     beta = rho / rho1
            #     # In-place update: bit-identical to 'p = z + beta * p' (same elementwise multiply
            #     # then add), but avoids one vector allocation + copy per iteration.
            #     # p = z + beta * p
            #     p *= beta
            #     p += z
            #
            # q = matvec(p)
            # denom = dot(xp, p, q)
            # denom_host = float(denom)
            #
            # if not bool(xp.isfinite(denom)) or denom_host <= 0.0:
            #     reason = "non-positive or nstop_atolon-finite p^T A p"
            #     break
            #
            # alpha = rho / denom
            #
            # if not bool(xp.isfinite(alpha)):
            #     reason = "non-finite alpha"
            #     break
            # Note on the blocking reads: the per-scalar reads and checks that used to be here are kept
            # below as comments. They required ~6 separate device->host transfers per iteration. Now
            # rho, denom and alpha are read back together with ONE transfer, and all decisions use those
            # host values in the same order. Nothing is written to x or r before these decisions, so the
            # breakdown exits return exactly the same x as before.
            if restart_p:
                p = z.copy() if hasattr(z, "copy") else z
                restart_p = False
            else:
                beta = rho / rho1
                # In-place update: bit-identical to 'p = z + beta * p' (same elementwise multiply
                # then add), but avoids one vector allocation + copy per iteration.
                # p = z + beta * p
                p *= beta
                p += z

            q = matvec(p)
            denom = dot(xp, p, q)
            alpha = rho / denom # device scalar, as before (keeps the dtype promotion of alpha * p unchanged)

            rho_host, denom_host, alpha_host = _host_scalars(xp, rho, denom, alpha)

            if not math.isfinite(rho_host):
                print("norm of r = ",norm(xp,r))
                print("norm of z=M@r= ",norm(xp,z))
                print("rho_host= ",rho_host)
                reason = "non-finite r^T M r"
                break

            if rho_host < 0.0:
                print("norm of r = ",norm(xp,r))
                print("norm of z=M@r= ",norm(xp,z))
                print("rho_host= ",rho_host)
                print("rho_host < 0 - going on as without preconditioner")
                reason = " r^T M r < 0"
                break

            #just proceed as nothing in the case it is exactly 0
            if rho_host == 0.0:
                print("rho_host == 0 - simply proceeding even though this should not happen")

            if not math.isfinite(denom_host) or denom_host <= 0.0:
                reason = "non-positive or nonfinite p^T A p"
                break

            if not math.isfinite(alpha_host):
                reason = "non-finite alpha"
                break

            # In-place updates: bit-identical to 'x = x + alpha*p' / 'r = r - alpha*q'
            # (same elementwise multiply then add/sub), but avoid two vector allocations
            # + copies per iteration (x/r are owned locally, safe to mutate).
            # x = x + alpha * p
            # r = r - alpha * q
            x += alpha * p
            r -= alpha * q

            iters += 1

            if callback is not None:
                callback(x)

            # Cheap recursive residual norm.
            resid = norm(xp, r)
            resid_host = float(resid)
            recursive_residual_norms.append(resid_host)

            # Original (kept for reference): if not bool(xp.isfinite(resid)):
            if not math.isfinite(resid_host): # same value as isfinite(resid), read back above
                reason = "non-finite recursive residual"
                break

            # Compute exact residual when the recursive residual claims convergence
            # or at the requested reliability-check frequency.
            need_exact_residual = resid_host <= stop_atol

            if recompute_residual_every:
                need_exact_residual = need_exact_residual or (iters % recompute_residual_every == 0)

            if need_exact_residual:
                r_exact = b - matvec(x)
                resid_exact = norm(xp, r_exact)
                resid_exact_host = float(resid_exact)
                true_residual_norms.append(resid_exact_host)

                # Original (kept for reference): if not bool(xp.isfinite(resid_exact)):
                if not math.isfinite(resid_exact_host): # same value as isfinite(resid_exact), read back above
                    reason = "non-finite true residual"
                    break

                if resid_exact_host <= stop_atol:
                    resid = resid_exact
                    resid_host = resid_exact_host
                    converged = True
                    reason = "true residual satisfies tolerance"
                    break

                gap = norm(xp, r_exact - r) / max(resid_exact_host, tiny)
                gap_host = float(gap)
                residual_gap_ratios.append(gap_host)

                if replace_residual and gap_host > residual_gap_tol:
                    r = r_exact
                    resid = resid_exact
                    resid_host = resid_exact_host
                    replacement_iterations.append(iters)

                    if restart_on_residual_replacement:
                        restart_p = True
                        rho = 0
                        continue

        info_int = 0 if converged else iters

        if return_info_dict:
            return x, {
                "converged": converged,
                "iterations": iters,
                "info": info_int,
                "reason": reason,
                "residual_norm": resid_host,
                "tolerance": stop_atol,
                "recursive_residual_norms": recursive_residual_norms,
                "true_residual_norms": true_residual_norms,
                "residual_gap_ratios": residual_gap_ratios,
                "replacement_iterations": replacement_iterations,
            }

        return x, info_int
