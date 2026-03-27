from numpy import double
import numpy as np
import KratosMultiphysics as KM

# Einsum optimization configuration
opt_type = "greedy"

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
        import numpy as np
        import scipy.sparse as sp
        import scipy.sparse.linalg as sp_linalg

        xp = np
        sparse = sp
        sparse_linalg = sp_linalg
        def _asnumpy(x): return np.asarray(x)
        asnumpy = _asnumpy
        USE_CUPY = False

    elif parallel_type == "GPU" or parallel_type == "gpu":
        import cupy as cp
        import cupyx.scipy.sparse as sp
        import cupyx.scipy.sparse.linalg as sp_linalg
        
        xp = cp
        sparse = sp
        sparse_linalg = sp_linalg
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

    def GetShapeFunctionsOnGaussPoints(self, dim: int, integration_order: int):
        if integration_order == 1:
            N = xp.array([[xp.ones(dim+1) / (dim+1)]])
        elif integration_order == 2:
            if dim == 2:
                N = xp.array([
                    [2.0/3.0, 1.0/6.0, 1.0/6.0],
                    [1.0/6.0, 2.0/3.0, 1.0/6.0],
                    [1.0/6.0, 1.0/6.0, 2.0/3.0]
                ])
            elif dim == 3:
                N = xp.array([
                    [0.5854101966249685, 0.1381966011250105, 0.1381966011250105, 0.1381966011250105],
                    [0.1381966011250105, 0.5854101966249685, 0.1381966011250105, 0.1381966011250105],
                    [0.1381966011250105, 0.1381966011250105, 0.5854101966249685, 0.1381966011250105],
                    [0.1381966011250105, 0.1381966011250105, 0.1381966011250105, 0.5854101966249685]
                ])
            else:
                raise Exception("Dimension not supported.")
        else:
            raise Exception("Integration order not supported.")

        return N

    def GetGaussIntegrationWeights(self, dim: int, integration_order: int):
        if integration_order == 1:
            if dim == 2:
                w = xp.array([1.0/2.0])
            elif dim == 3:
                w = xp.array([1.0/6.0])
            else:
                raise Exception("Dimension not supported.")
        elif integration_order == 2:
            if dim == 2:
                w = xp.array([
                    1.0/6.0,
                    1.0/6.0,
                    1.0/6.0
                ])
            elif dim == 3:
                w = xp.array([
                    1.0/24.0,
                    1.0/24.0,
                    1.0/24.0,
                    1.0/24.0
                ])
            else:
                raise Exception("Dimension not supported.")
        else:
            raise Exception("Integration order not supported.")

        return w

    def ComputeElementalDivergence(self, DN: np.ndarray, uel: np.ndarray):
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
        """
        return xp.einsum("nij,nij->n", DN, uel, optimize=opt_type)

    def ComputeElementwiseNodalDivergence(self, N: np.ndarray, DN: np.ndarray, uel: np.ndarray):
        """
        Computes the nodal weighted divergence of a vector field u.

        This term is (w, ∇·u).
        Using Einstein notation: out[e,i] = sum_k sum_l N[i]*DN[e,k,l]*uel[e,k,l]

        Parameters
        ----------
        N : ndarray
            Numpy array with shape (n_in_el,).
        DN : ndarray
            Numpy array with shape (Nelem, n_in_el, dim).
        uel : ndarray
            Numpy array with shape (Nelem, n_in_el, dim).
        out : ndarray
            Output array, expected to have shape (Nelem, n_in_el).
        """
        n_in_el = DN.shape[1]
        if N.shape[0] != n_in_el or N.ndim!=1:
            raise Exception("wrong size of N")
        if uel.shape != DN.shape:
            raise Exception("wrong size of uel")

        return xp.einsum("I,eJk,eJk->eI",N,DN,uel,optimize=opt_type)

    def Compute_N_DN(self, N: np.ndarray, DN: np.ndarray, pel: np.ndarray):
        """
        Computes the term (w, ∇p).

        Using Einstein notation: out[e,I,k] = N[I]*DN[e,J,k]*pel[e,J]

        Parameters
        ----------
        N : ndarray
            Numpy array with shape (n_in_el,).
        DN : ndarray
            Numpy array with shape (Nelem, n_in_el, dim).
        pel : ndarray
            Numpy array with shape (Nelem, n_in_el), i.e., one scalar per node of the element.
        """
        nelem = DN.shape[0]
        n_in_el = DN.shape[1]

        if pel.shape != (nelem, n_in_el):
            raise ValueError("pel must have shape (nelem, n_in_el) for scalar case. Current shape is:",field.shape)

        return xp.einsum("I,eJk,eJ->eIk", N, DN, pel, optimize=opt_type)

    def Compute_DN_N(self, N: np.ndarray, DN: np.ndarray, pel: np.ndarray):
        """
        Computes the term (∇·w, p).

        Using Einstein notation: out[e,I,k] = sum_J DN[e,I,k]*pel[e,J]*N[J]

        Parameters
        ----------
        N : ndarray
            Numpy array with shape (n_in_el,).
        DN : ndarray
            Numpy array with shape (Nelem, n_in_el, dim).
        pel : ndarray
            Numpy array with shape (Nelem, n_in_el), i.e., one scalar per node of the element.
        """
        return xp.einsum("eIk,J,eJ->eIk", DN, N, pel, optimize=opt_type)

    def ComputeLaplacianMatrix(self, DN: np.ndarray):
        """
        Computes the Laplacian local matrix term (∇N, ∇N).

        Using Einstein notation: out[e,i,j] = DN[e,i,m]*DN[e,j,m]

        Parameters
        ----------
        DN : ndarray
            Numpy array with shape (Nelem, n_in_el, dim).
        out : ndarray
            Output array, expected to have shape (Nelem, n_in_el, n_in_el).
        """
        return xp.einsum("eim,ejm->eij", DN, DN, optimize=opt_type)

    def ApplyLaplacian(self, DN: np.ndarray, field: np.ndarray):
        """
        Computes the term (∇q, ∇field) - scalar field
                       or (∇w, ∇field) - vector field

        Using Einstein notation:
        if field is scalar: out[e,i] = DN[e,i,m]*DN[e,j,m]*field[e,j]
        if field is vector: out[e,i,k] = DN[e,i,m]*DN[e,j,m]*field[e,j,k]

        Parameters
        ----------
        DN : ndarray
            Numpy array with shape (Nelem, n_in_el, dim).
        field : ndarray
            Numpy array with shape (Nelem, n_in_el). - scalar case
            or
            Numpy array with shape (Nelem, n_in_el, dim). - vector case
        out : ndarray
            Output array, expected to have shape (Nelem, n_in_el).
            or
            Output array, expected to have shape (Nelem, n_in_el, dim).
        """

        if field.ndim == 2: # scalar case
            return xp.einsum("eIm,eJm,eJ->eI", DN, DN, field, optimize=opt_type)
        elif field.ndim == 3: # vector case
            return xp.einsum("eIm,eJm,eJk->eIk", DN, DN, field, optimize=opt_type)

        raise ValueError("Field must have shape (nelem,nnode) or (nelem,nnode,dim)")

    def InterpolateValue(self, N, field):
        """
        Interpolates field values at Gauss points.

        Parameters
        ----------
        N : (nnode,) or (ngauss, nnode)
        field : (nelem, nnode) or (nelem, nnode, dim)

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

        # Scalar field
        if field.ndim == 2:
            # (ngauss, nnode) x (nelem, nnode) -> (nelem, ngauss)
            return xp.einsum("gn,en->eg", N, field, optimize=opt_type)

        # Vector field
        if field.ndim == 3:
            # (ngauss, nnode) x (nelem, nnode, dim) -> (nelem, ngauss, dim)
            return xp.einsum("gn,end->egd", N, field, optimize=opt_type)

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

    def ComputeElementalGradient(self, DN: np.ndarray, field: np.ndarray):
        """
        Compute the gradient of a scalar or vector field using element-dependent DN.
        that is for every element we compute:
            grad[k]   = sum_I DN_I/Dx_k p_I  - scalar case
            grad[k,l] = sum_I DN_I/Dx_l v_Ik - vector case
        Parameters
        ----------
        DN : (nelem, nnode, ndim)
            Shape function derivatives for each element.

        field : array
            Field values at element nodes. Must have shape:
            - (nelem, nnode)          for a scalar field
            - (nelem, nnode, ncomp)   for a vector field

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

            return xp.einsum('ei,eik->ek', field, DN, optimize=opt_type)

        # ------------------------------
        # Vector field
        # ------------------------------
        elif field.ndim == 3:
            if field.shape[0] != nelem or field.shape[1] != nnode:
                raise ValueError("field must have shape (nelem, nnode, ncomp). Current shape is:",field.shape)

            return xp.einsum('eik,eil->ekl', field, DN, optimize=opt_type)

        # ------------------------------
        # Invalid field
        # ------------------------------
        raise ValueError("field must have 2 dims (scalar) or 3 dims (vector), Current shape of field is:",field.shape)

    def ComputeBodyForceContribution(self, N, b_elemental):
        """
        Compute (w, b) contribution for multiple Gauss points.

        Parameters
        ----------
        N : (ngauss, nnode)
            shape function values at Gauss points

        b : (nelem, nnode, dim)
            body force at Gauss points

        Returns
        -------
        (nelem, ngauss, dim)
        """

        return xp.einsum("gn,end->egd", N, b_elemental, optimize=opt_type)

    def ComputeConvectiveContribution(self, N, grad_u, a_gauss):
        """
        Compute (w,a·∇u) although with the definition we employ for ∇u this is actually (w,∇u·a)
        Note that this function assumes ∇u to be constant within the element

        Parameters
        ----------
        N : (ngauss, nnode)
        grad_u : (nelem, ndim, ndim) or (nelem, ndim)
        a_gauss : (nelem, ngauss, ndim)

        Returns
        -------
        vector case: (nelem, ngauss, nnode, ndim) and (nelem, ngauss, ndim)
        scalar case: (nelem, ngauss, nnode) and (nelem, ngauss, ndim)
        """

        if grad_u.ndim == 3:
            
            conv = xp.einsum('ekl,egl->egk', grad_u, a_gauss, optimize=opt_type) #a·∇u
            out = N[None, :, :, None] * conv[:, :, None, :] #w,a·∇u

        elif grad_u.ndim == 2:
            conv = xp.einsum('el,egl->eg', grad_u, a_gauss, optimize=opt_type) #a·∇u
            out = N[None, :, :] * conv[:, :, None] #w,a·∇u

        else:
            raise ValueError("grad_u must have 2 dims (scalar) or 3 dims (vector)")

        return out

    def ComputeMomentumStabilization(self, N, DN, u_elemental, a_gauss, pi_gauss, rho):
        """
        Compute convection + convective stabilization in a single pass.

        (w, a·∇u) + (ρ a·∇w, ρ a·∇u) - (ρ a·∇w, Π)

        Parameters
        ----------
        N : (ngauss, nnode)
        DN : (nelem, nnode, ndim)
        u_elemental : (nelem, nnode, ndim)
        a_gauss : (nelem, ngauss, ndim)
        pi_gauss : (nelem, ngauss, ndim)
        rho : float

        Returns
        -------
        out : (nelem, ngauss, nnode, ndim)
            Total convective contribution (Galerkin + stabilization)
        """

        # --- discrete operator: a · ∇N ---
        adv = xp.einsum("egd,end->egn", a_gauss, DN, optimize=opt_type)  # (E,G,nnode)

        # --- a·∇u = Σ_J adv_J u_J ---
        conv = xp.einsum("egn,end->egd", adv, u_elemental, optimize=opt_type)  # (E,G,dim)

        # --- Galerkin term ---
        out = N[None, :, :, None] * conv[:, :, None, :]  # (E,G,nnode,dim)

        # --- stabilization ---
        beta = rho * adv                                # (E,G,nnode)

        tmp = xp.einsum("egn,end->egd", beta, u_elemental, optimize=opt_type)

        out += beta[:, :, :, None] * (tmp[:, :, None, :] - pi_gauss[:, :, None, :])

        return out

    def ComputeDivDivStabilization(self, N: np.array, DN: np.ndarray, u_elemental : np.ndarray, Pi_div_elemental: np.ndarray):
        output = xp.einsum("eik,ejl,ejl->eik", DN, DN, u_elemental, optimize=opt_type)
        elem_scratch = xp.einsum("eik,j,ej->eik", DN, N, Pi_div_elemental, optimize=opt_type)
        output -= elem_scratch
        return output

    def ComputePressureStabilization_ProjectionTerm(self, N: np.ndarray, DN: np.ndarray, Pi_press_el: np.ndarray):
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
