import cupy
#import pyamg
import numpy as np
import scipy.sparse as sp
import KratosMultiphysics as KM
import KratosMultiphysics.FluidDynamicsApplication as CFDApp

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

class CFDUtils:

    def __init__(self):
        if not _configured:
            raise RuntimeError("Backend not set. Call 'configure()' first.")
        self._tmp_nelem = 0
        self._tmp_n_in_el = 0
        self._tmp_dim = 0

    def allocate_temporaries(self, nelem, n_in_el, dim):
        """
        Pre-allocate shared temporary buffers used by compute methods.

        These buffers are reused across calls under the assumption that
        only one compute method executes at a time.  The returned arrays
        from compute methods are views into these buffers and remain
        valid only until the next call.

        Parameters
        ----------
        nelem : int
            Number of elements.
        n_in_el : int
            Number of nodes per element.
        dim : int
            Spatial dimension.
        """
        self._tmp_nelem = nelem
        self._tmp_n_in_el = n_in_el
        self._tmp_dim = dim

        nd = n_in_el * dim

        # Scalar per element — shape (nelem,)
        self._tmp_e = xp.empty(nelem, dtype=PRECISION)

        # Flat node×dim per element — shape (nelem, n_in_el*dim)
        self._tmp_e_nd = xp.empty((nelem, nd), dtype=PRECISION)

        # Per-element dim vector — shape (nelem, dim)
        self._tmp_ed = xp.empty((nelem, dim), dtype=PRECISION)

        # Per-element node×dim — shape (nelem, n_in_el, dim)
        self._tmp_eik = xp.empty((nelem, n_in_el, dim), dtype=PRECISION)

        # Per-element node scalar — shape (nelem, n_in_el)
        self._tmp_en = xp.empty((nelem, n_in_el), dtype=PRECISION)

        # Per-element dim×dim — shape (nelem, dim, dim)
        self._tmp_edd = xp.empty((nelem, dim, dim), dtype=PRECISION)

    def _ensure_temporaries(self, nelem, n_in_el, dim):
        """Lazily allocate temporaries on first use or if dimensions change."""
        if nelem != self._tmp_nelem or n_in_el != self._tmp_n_in_el or dim != self._tmp_dim:
            self.allocate_temporaries(nelem, n_in_el, dim)

    def GetShapeFunctionsOnGaussPoints(self, dim: int, integration_order: int):
        if integration_order == 1:
            N = xp.array([[xp.ones(dim+1) / (dim+1)]])
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

        # Element-wise product into pre-allocated buffer
        xp.multiply(DN_flat, uel_flat, out=self._tmp_e_nd)

        # Sum rows into caller-provided (or freshly allocated) output
        xp.sum(self._tmp_e_nd, axis=1, out=out)

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
        self._ensure_temporaries(nelem, n_in_el, dim)
        nd = n_in_el * dim

        if out is None:
            out = xp.empty((nelem, n_in_el), dtype=DN.dtype)

        # Step 1: divergence div[e] = sum_{J,k} DN[e,J,k] * uel[e,J,k]
        DN_flat = DN.reshape(nelem, nd)
        uel_flat = uel.reshape(nelem, nd)
        xp.multiply(DN_flat, uel_flat, out=self._tmp_e_nd)
        xp.sum(self._tmp_e_nd, axis=1, out=self._tmp_e)

        # Step 2: out[e,I] = N[I] * div[e]
        xp.multiply(self._tmp_e[:, None], N[None, :], out=out)

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
            xp.matmul(a_elemental, grad_u.swapaxes(1, 2), out=out)
            return out

        elif grad_u.ndim == 2: # scalar field
            if out is None:
                out = xp.empty((nelem, n_in_el), dtype=PRECISION)

            # (E, I, L) @ (E, L, 1) -> (E, I, 1) -> squeeze to (E, I)
            xp.matmul(a_elemental, grad_u[:, :, None], out=out[:, :, None])
            return out

        else:
            raise ValueError("grad_u must have 2 dims (scalar) or 3 dims (vector)")

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
            raise ValueError("pel must have shape (nelem, n_in_el) for scalar case. Current shape is:",field.shape)

        self._ensure_temporaries(nelem, n_in_el, dim)

        if out is None:
            out = xp.empty((nelem, n_in_el, dim), dtype=DN.dtype)

        # Step 1: grad_p[e,k] = sum_J DN[e,J,k] * pel[e,J]
        xp.matmul(DN.swapaxes(1, 2), pel[:, :, None], out=self._tmp_ed[:, :, None])

        # Step 2: out[e,I,k] = N[I] * grad_p[e,k]
        xp.multiply(N[None, :, None], self._tmp_ed[:, None, :], out=out)

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
        xp.dot(pel, N, out=self._tmp_e)

        # Step 2: out[e,I,k] = DN[e,I,k] * s[e]
        xp.multiply(DN, self._tmp_e[:, None, None], out=out)

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
        self._ensure_temporaries(nelem, nnode, dim)

        # Scalar field
        if field.ndim == 2:
            if out is None:
                out = xp.empty((nelem, nnode), dtype=field.dtype)

            # Use _tmp_eik for intermediate (shape: nelem, nnode, dim)
            # grad[:, m] = sum_node DN[:, :, m] * field
            xp.multiply(DN, field[:, :, None], out=self._tmp_eik)
            xp.sum(self._tmp_eik, axis=1, out=self._tmp_ed)  # (nelem, dim)

            # Accumulate: out[e, i] = sum_m DN[e, i, m] * grad[e, m]
            # This is: out = DN @ grad[:, :, None] → (nelem, nnode, 1)
            xp.matmul(DN, self._tmp_ed[:, :, None], out=out[:, :, None])

            return out

        # Vector field
        elif field.ndim == 3:
            _, _, ncomp = field.shape
            if out is None:
                out = xp.empty((nelem, nnode, ncomp), dtype=field.dtype)

            # Use _tmp_edd for grad (shape: nelem, dim, dim) — slice to (nelem, dim, ncomp)
            grad_view = self._tmp_edd[:, :, :ncomp]

            # grad = DN^T @ field → (nelem, dim, ncomp)
            xp.matmul(DN.swapaxes(1, 2), field, out=grad_view)

            # out = DN @ grad → (nelem, nnode, ncomp)
            xp.matmul(DN, grad_view, out=out)

            return out

        raise ValueError("Field must have shape (E,nnode) or (E,nnode,dim)")

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
            xp.matmul(field.swapaxes(1, 2), DN, out=out) #TODO: verify if this is efficient
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
        xp.matmul(M_e, b_elemental, out=out)   # (e, node, dim)

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

            xp.matmul(M_e, elem_conv, out=out)
            return out

        elif elem_conv.ndim == 2: # scalar field
            nelem, n_node = elem_conv.shape
            dim = n_node - 1
            M_e = self.GetElementalMassMatrix(dim) #elemental mass matrix (n_node, n_node)

            if out is None:
                out = xp.empty((nelem, n_node), dtype=elem_conv.dtype)
            elif out.shape != (nelem, n_node):
                raise ValueError("out must have shape (n_elem, n_node)")

            xp.matmul(elem_conv, M_e.T, out=out)
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

        beta = conv_elemental - pi_elemental
        n_elem, n_node, n_dim = beta.shape

        if out is None:
            out = xp.empty((n_elem, n_node, n_dim), dtype=beta.dtype)
        elif out.shape != (n_elem, n_node, n_dim):
            raise ValueError(
                f"out must have shape (n_elem, n_node, n_dim) = ({n_elem}, {n_node}, {n_dim}), "
                f"but got shape {out.shape}"
            )
        # Note: out is fully overwritten by xp.matmul(DN, S, out=out) on line 828,
        # so no zeroing is needed when reusing pre-allocated arrays.

        M_e = self.GetElementalMassMatrix(n_dim)
        
        # Original computation (commented for reference):
        # S = xp.zeros((n_elem, n_dim, n_dim), dtype=beta.dtype)
        # for j in range(n_node):
        #     for m in range(n_node):
        #         w = M_e[j, m]
        #         S += w * (a_elemental[:, j, :, None] * beta[:, m, None, :])
        # out = xp.matmul(DN, S, out=out)
        
        # Vectorized computation (no loops, no 4th-order tensor):
        # out[e,i,l] = sum_{j,m,k} DN[e,i,k] * M[j,m] * a[e,j,k] * beta[e,m,l]
        # Factor as: S[e,k,l] = sum_{j,m} M[j,m] * a[e,j,k] * beta[e,m,l]
        # Then: out[e,i,l] = sum_k DN[e,i,k] * S[e,k,l]
        
        self._ensure_temporaries(n_elem, n_node, n_dim)
        
        # Step 1: Compute temp[e,j,l] = sum_m M[j,m] * beta[e,m,l]
        # Batched matmul: M @ beta[e] for all elements at once
        # Use _tmp_eik as temp buffer (shape: nelem, n_node, n_dim)
        xp.matmul(M_e[None, :, :], beta, out=self._tmp_eik)
        
        # Step 2: Compute S[e,k,l] = sum_j a[e,j,k] * temp[e,j,l]
        # Transpose a to (n_elem, n_dim, n_node) then matmul with temp
        # Use _tmp_edd as S buffer (shape: nelem, dim, dim)
        xp.matmul(a_elemental.swapaxes(1, 2), self._tmp_eik, out=self._tmp_edd)
        
        # Step 3: Compute out[e,i,l] = sum_k DN[e,i,k] * S[e,k,l]
        # Use out parameter to avoid implicit temporary from matmul
        xp.matmul(DN, self._tmp_edd, out=out)

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

        # Step 1: divergence of u → (E,)
        # Use _tmp_eik as intermediate buffer for DN * u_elemental
        xp.multiply(DN, u_elemental, out=self._tmp_eik)
        xp.sum(self._tmp_eik, axis=(1, 2), out=self._tmp_e)

        # Step 2: projection of Π → (E,)
        # Use _tmp_en as intermediate buffer for N[None, :] * Pi_div_elemental
        xp.multiply(N[None, :], Pi_div_elemental, out=self._tmp_en)
        xp.sum(self._tmp_en, axis=1, out=self._tmp_ed[:, 0])  # reuse _tmp_ed[:, 0] as (E,) view

        # Step 3: scalar residual → (E,) — in-place subtraction
        xp.subtract(self._tmp_e, self._tmp_ed[:, 0], out=self._tmp_e)

        # Step 4: apply DN → (E, nnode, dim)
        xp.multiply(DN, self._tmp_e[:, None, None], out=out)

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
        self._ensure_temporaries(nelem, n_node, n_dim)

        if out is None:
            out = xp.empty((nelem, n_node), dtype=DN.dtype)

        # Step 1: Compute temp[e,k] = sum_J N[J] * Pi_press_el[e,J,k]
        # Using matmul: Pi_press_el.swapaxes(1,2) @ N gives (nelem, dim, n_node) @ (n_node,) -> (nelem, dim)
        xp.matmul(Pi_press_el.swapaxes(1, 2), N, out=self._tmp_ed)

        # Step 2: Compute out[e,I] = sum_k DN[e,I,k] * temp[e,k]
        # Use _tmp_eik as intermediate buffer for DN * temp, then sum over k
        xp.multiply(DN, self._tmp_ed[:, None, :], out=self._tmp_eik)
        xp.sum(self._tmp_eik, axis=2, out=out)

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

