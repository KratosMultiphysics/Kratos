import pyamg
import numpy as np
import scipy.sparse as sp
import KratosMultiphysics as KM

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

    def ApplyLaplacian(self, DN: np.ndarray, field: np.ndarray):
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

        Returns
        -------
        (E, nnode) or (E, nnode, dim)
            Laplacian applied to the field.
        """

        nelem, _, dim = DN.shape

        # Scalar field
        if field.ndim == 2:
            grad = xp.zeros((nelem, dim), dtype=field.dtype)
            for m in range(dim):
                grad[:, m] = xp.sum(DN[:, :, m] * field, axis=1)

            out = xp.zeros_like(field)
            for m in range(dim):
                out += DN[:, :, m] * grad[:, m][:, None]

            return out

        # Vector field
        elif field.ndim == 3:
            _, _, ncomp = field.shape

            grad = xp.zeros((nelem, dim, ncomp), dtype=field.dtype)
            for m in range(dim):
                grad[:, m, :] = xp.sum(
                    DN[:, :, m][:, :, None] * field,
                    axis=1
                )

            out = xp.zeros_like(field)
            for m in range(dim):
                out += DN[:, :, m][:, :, None] * grad[:, m, :][:, None, :]

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

    def ComputeConvectiveContribution(self, elem_conv):
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
            _, _, n_dim = elem_conv.shape
            M_e = self.GetElementalMassMatrix(n_dim) #elemental mass matrix (n_node, n_node)
            return xp.einsum('ij,ejk->eik', M_e, elem_conv) # (n_elem, n_node, n_dim)

        elif elem_conv.ndim == 2: # scalar field
            _, n_dim = elem_conv.shape
            M_e = self.GetElementalMassMatrix(n_dim) #elemental mass matrix (n_node, n_node)
            return xp.einsum('ij,ej->ei', M_e, elem_conv) # (n_elem, n_node)

        else:
            raise ValueError("Provided elemental nodal convective operator must have 2 dims (scalar) or 3 dims (vector)")

        return _

    def ComputeMomentumStabilization(self, DN, a_elemental, conv_elemental, pi_elemental):
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

        Returns
        -------
        (n_elem, n_node, n_dim)

        Notes
        -----
        - Equivalent to Gauss integration under consistent interpolation.
        - Geometric scaling (detJ / volume) applied externally.
        """

        beta = conv_elemental - pi_elemental
        n_elem, n_node, n_dim = beta.shape

        # --- build S[k,l] ---
        M_e = self.GetElementalMassMatrix(beta.shape[2])
        S = xp.zeros((n_elem, n_dim, n_dim), dtype=beta.dtype)
        for j in range(n_node):
            for m in range(n_node):
                w = M_e[j, m]
                S += w * (a_elemental[:, j, :, None] * beta[:, m, None, :]) # a_j,k * beta_m,l

        # --- contract with DN ---
        R = xp.zeros((n_elem, n_node, n_dim), dtype=beta.dtype)
        for k in range(n_dim):
            R += DN[:, :, k][:, :, None] * S[:, k, :][:, None, :]

        return R

    def ComputeDivDivStabilization(self, N: np.array, DN: np.ndarray, u_elemental : np.ndarray, Pi_div_elemental: np.ndarray):
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

        Returns
        -------
        (E, nnode, dim)
        """

        # Step 1: divergence of u → (E,)
        div_u = xp.sum(DN * u_elemental, axis=(1, 2))

        # Step 2: projection of Π → (E,)
        proj = xp.sum(N[None, :] * Pi_div_elemental, axis=1)

        # Step 3: scalar residual → (E,)
        div_res = div_u - proj

        # Step 4: apply DN → (E, nnode, dim)
        return DN * div_res[:, None, None]

    def ComputePressureStabilizationProjectionTerm(self, N: np.ndarray, DN: np.ndarray, Pi_press_el: np.ndarray):
        """
        implements (∇q,Pi_pressure)

        this corresponds in einstain notation to
        out[i,k] = sum_l N_I ∇u_kl a_l
        on every element

        Parameters
        ----------
        N : (nelem, n_in_el)
            shape function values at the gauss point

        grad_u : (nelem, ndim, ndim) in the vector case or (nelem, ndim) in the scalar case
            gradient of the velocity at the gauss point

        a: (nelem,ndim)
            convective velocity on the gauss point

        out: (nelem, n_in_el, ndim)

        """
        return xp.einsum("eIk,J,eJk->eI",DN,N,Pi_press_el, optimize=opt_type)

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

class GraphBasedSA_AMG:
    def __init__(self, A_pattern_gpu, max_coarse=300):
        self.shape = A_pattern_gpu.shape
        self.dtype = A_pattern_gpu.dtype
        self.max_coarse = max_coarse

        # ---------------------------------------------------------
        # CPU PHASE: Purely structural graph analysis
        # ---------------------------------------------------------
        indptr_cpu = A_pattern_gpu.indptr.get()
        indices_cpu = A_pattern_gpu.indices.get()

        dummy_data = A_pattern_gpu.data.get()
        A_dummy_cpu = sp.csr_matrix((dummy_data, indices_cpu, indptr_cpu), shape=self.shape)

        ml_dummy = pyamg.smoothed_aggregation_solver(
            A_dummy_cpu,
            max_coarse=max_coarse,
            strength=('symmetric', {'theta': 0.0}),
            keep=True
        )

        # ---------------------------------------------------------
        # GPU TRANSFER PHASE
        # ---------------------------------------------------------
        self.T_operators_gpu = []
        for i in range(len(ml_dummy.levels) - 1):
            T_cpu = ml_dummy.levels[i].T.tocsr()
            # Send T permanently to the GPU
            T_gpu = sparse.csr_matrix(T_cpu)
            self.T_operators_gpu.append(T_gpu)

        self.levels = []

    def update_matrix_values(self, A_new):
        """
        Phase 2: GPU only - Fast rebuild of P, R, and A_coarse.
        """
        self.levels = []
        A_current = A_new #.tocsr()

        # Lock omegas to prevent Float64 upcasting
        #omega_prolong = xp.float32(4.0 / 3.0) if self.dtype == xp.float32 else 4.0 / 3.0
        omega_prolong = PRECISION(4.0/3.0)

        for lvl_idx, T in enumerate(self.T_operators_gpu):
            lvl_dict = {'A': A_current}

            # Let CuPy find and extract the diagonal natively on the GPU
            diag = A_current.diagonal()

            # Protect against near-zeros causing float32 explosions
            #diag[xp.abs(diag) < 1e-7] = 1.0
            invD = xp.reciprocal(diag) #
            #invD = PRECISION(1.0) / diag #TODO: check precision!
            lvl_dict['invD'] = invD

            # Build Smoothed Prolongator using @
            A_T = A_current @ T
            scaled_A_T = sparse.spdiags(invD, 0, A_current.shape[0], A_current.shape[1]) * A_T

            P = T - (scaled_A_T * omega_prolong)
            lvl_dict['P'] = P
            R = P.T
            self.levels.append(lvl_dict)

            # Generate the new structure natively on GPU using @
            A_current = R @ A_current @ P

        # --- Base (Coarsest) Level Handling ---
        base_dict = {'A': A_current}

        if A_current.shape[0] <= self.max_coarse:
            base_dict['is_dense_solve'] = True
            A_dense_64 = A_current.todense().astype(xp.float64)
            # Use pseudo-inverse to survive singular coarse grids
            A_inv_64 = xp.linalg.pinv(A_dense_64)
            base_dict['A_dense_inv'] = A_inv_64.astype(self.dtype)

        else:
            base_dict['is_dense_solve'] = False

            # Native extraction for the base level
            diag = A_current.diagonal()
            diag[xp.abs(diag) < 1e-7] = 1.0
            base_dict['invD'] = PRECISION(1.0) / diag

        self.levels.append(base_dict)

    def _v_cycle(self, b, level_idx=0):
        lvl = self.levels[level_idx]
        A = lvl['A']

        # Base Case Solve
        if level_idx == len(self.levels) - 1:
            if lvl['is_dense_solve']:
                return lvl['A_dense_inv'] @ b
            else:
                x = xp.zeros_like(b)
                omega_jacobi = PRECISION(0.67)
                for _ in range(20):
                    r = b - (A @ x)
                    r *= lvl['invD']
                    r *= omega_jacobi
                    x += r
                return x

        P = lvl['P']
        invD = lvl['invD']
        omega = PRECISION(0.67)
        sweeps = 2

        x = xp.zeros_like(b)

        # Pre-smoothing
        for _ in range(sweeps):
            r = b - (A @ x)
            r *= invD
            r *= omega
            x += r

        # Restrict
        r = b - (A @ x)
        b_coarse = P.T @ r

        # Coarse grid correction
        e_coarse = self._v_cycle(b_coarse, level_idx + 1)

        # Prolongate
        x += P @ e_coarse

        # Post-smoothing
        for _ in range(sweeps):
            r = b - (A @ x)
            r *= invD
            r *= omega
            x += r

        return x

    def matvec(self, b):
        return self._v_cycle(b, level_idx=0)

    def aspreconditioner(self):
        return sparse_linalg.LinearOperator(shape=self.shape, matvec=self.matvec, dtype=self.dtype)

