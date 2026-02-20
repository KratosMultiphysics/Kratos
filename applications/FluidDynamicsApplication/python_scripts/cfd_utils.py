from numpy import double
import numpy as np
import KratosMultiphysics as KM

# Configuration flag for backend selection
# This can be changed to True to use CuPy if available
USE_CUPY = True

# Precision configuration
PRECISION = np.float64

if USE_CUPY:
    try:
        import cupy as xp
        import cupyx.scipy.sparse as sparse
        import cupyx.scipy.sparse.linalg as sparse_linalg
        asnumpy = xp.asnumpy
    except ImportError:
        print("CuPy not found. Falling back to NumPy.")
        import numpy as xp
        import scipy.sparse as sparse
        import scipy.sparse.linalg as sparse_linalg
        def asnumpy(x):
            return np.asarray(x)
        USE_CUPY = False
else:
    import numpy as xp
    import scipy.sparse as sparse
    import scipy.sparse.linalg as sparse_linalg
    def asnumpy(x):
        return np.asarray(x)

class CFDUtils:
    def __init__(self):
        #allocate auxiliary arrays
        self.aux_array1 = xp.empty(0, dtype=PRECISION)

    def ComputeElementalDivergence(self, DN, uel, out):
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
        out : ndarray
            Output array, expected to have shape (Nelem,).
        """
        nelem = DN.shape[0]

        #verify shape of output array
        if(out.shape != (nelem,)):
            raise ValueError("Output array has wrong shape")

        #store element by element divergence in aux_array1
        xp.einsum("nij,nij->n", DN, uel, out=out, optimize=True)

    def ComputeElementwiseNodalDivergence(self, N, DN, uel, out):
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
        nelem = DN.shape[0]
        n_in_el = DN.shape[1]
        dim = DN.shape[2]
        if N.shape[0] != n_in_el or N.ndim!=1:
            raise Exception("wrong size of N")
        if uel.shape != DN.shape:
            raise Exception("wrong size of uel")

        xp.einsum("I,eJk,eJk->eI",N,DN,uel,out=out,optimize=True)

    def Compute_N_DN(self, N, DN, pel, out):
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
        out : ndarray
            Output array, expected to have shape (Nelem, n_in_el).
        """
        nelem = DN.shape[0]
        n_in_el = DN.shape[1]
        ndim = DN.shape[2]

        if pel.shape != (nelem, n_in_el):
            raise ValueError("pel must have shape (nelem, n_in_el) for scalar case. Current shape is:",field.shape)

        if out.shape != (nelem, n_in_el,ndim):
            raise ValueError(f"out must have shape (nelem, {ndim}) for scalar case. Current shape is:",out.shape)

        xp.einsum("I,eJk,eJ->eIk", N,DN, pel,  out=out, optimize=True)

    def Compute_DN_N(self, N, DN, pel, out):
        """
        Computes the term (∇w, p).

        Using Einstein notation: out[e,I,k] = sum_J DN[e,I,k]*pel[e,J]*N[J]

        Parameters
        ----------
        N : ndarray
            Numpy array with shape (n_in_el,).
        DN : ndarray
            Numpy array with shape (Nelem, n_in_el, dim).
        pel : ndarray
            Numpy array with shape (Nelem, n_in_el), i.e., one scalar per node of the element.
        out : ndarray
            Output array, expected to have shape (Nelem, n_in_el, dim).
        """
        nelem = out.shape[0]

        xp.einsum("eIk,J,eJ->eIk",DN,N,pel,out=out)

        # # aux_array1 must have shape (nelem,)
        # if self.aux_array1.shape != (nelem,):
        #     self.aux_array1 = np.empty(nelem, dtype=PRECISION)

        # # 1) aux[e] = pel[e,k] * N[k]
        # np.einsum("ek,k->e", pel, N, out=self.aux_array1)

        # # 2) out[e,i,m] = DN[e,i,m] * aux[e]
        # out[:] = DN
        # out *= self.aux_array1[:, np.newaxis, np.newaxis]

    def ComputeLaplacianMatrix(self, DN, out):
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
        xp.einsum("eim,ejm->eij", DN, DN, out=out, optimize=True)

    def ApplyLaplacian(self, DN, field, out):
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

            xp.einsum("eIm,eJm,eJ->eI", DN, DN, field, out=out)

        elif field.ndim == 3:

            xp.einsum("eIm,eJm,eJk->eIk", DN, DN, field, out=out)


    def InterpolateValue(self, N, field, out):
        """
        Computes:
            pgauss[e]     = sum_k N[k] * pel[e,k]
            vgauss[e,i]   = sum_k N[k] * vel[e,k,i]
        depending on the shape of 'field'.

        Parameters
        ----------
        N : (nnode,)
            Shape function values at the Gauss point.
        field : (nelem, nnode) or (nelem, nnode, dim)
            Element field to project (pel or vel).
        out : preallocated array
            Must have shape (nelem,) or (nelem, dim) accordingly.
        """

        nelem = field.shape[0]
        nnode = N.shape[0]

        # Scalar case: pel[e,k] → pgauss[e]
        if field.ndim == 2:
            if out.shape != (nelem,):
                raise ValueError(f"Output must have shape ({nelem},)")

            # pgauss[e] = N[k] * pel[e,k]
            xp.einsum("k,ek->e", N, field, out=out)
            return

        # Vector case: vel[e,k,i] → vgauss[e,i]
        if field.ndim == 3:
            dim = field.shape[2]
            if out.shape != (nelem, dim):
                raise ValueError(f"Output must have shape ({nelem},{dim})")

            # vgauss[e,i] = N[k] * vel[e,k,i]
            xp.einsum("k,eki->ei", N, field, out=out)
            return

        raise ValueError("Field must have shape (nelem,nnode) or (nelem,nnode,dim)")

    def AssembleVector(self, conn, vals, out):
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
        ###NOOOOOOOOOOOOOOO
        xp.add.at(out, conn, vals)

    def ComputeElementalGradient(self, DN, field, out):
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

        out : array
            Preallocated output buffer. Must have shape:
            - (nelem, ndim)           for scalar field
            - (nelem, ncomp, ndim)    for vector field

        Notes
        -----
        Computes:
            Scalar: out[e,k]      = sum_i field[e,i]      * DN[e,i,k]
            Vector: out[e,j,k]    = sum_i field[e,i,j]    * DN[e,i,k]

        No temporaries are allocated; `out` is written in-place.
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

            if out.shape != (nelem, ndim):
                raise ValueError(f"out must have shape (nelem, {ndim}) for scalar case. Current shape is:",out.shape)

            # grad[e,k] = sum_i field[e,i] * DN[e,i,k]
            xp.einsum('ei,eik->ek', field, DN, out=out)
            return

        # ------------------------------
        # Vector field
        # ------------------------------
        elif field.ndim == 3:
            if field.shape[0] != nelem or field.shape[1] != nnode:
                raise ValueError("field must have shape (nelem, nnode, ncomp). Current shape is:",field.shape)

            ncomp = field.shape[2]

            if out.shape != (nelem, ncomp, ndim):
                raise ValueError(
                    f"out must have shape (nelem, {ncomp}, {ndim}) for vector case. Current shape is:",field.shape
                )

            # grad[e,k,l] = sum_i field[e,i,k] * DN[e,i,l]
            xp.einsum('eik,eil->ekl', field, DN, out=out)
            return

        # ------------------------------
        # Invalid field
        # ------------------------------
        print(field.ndim)
        raise ValueError("field must have 2 dims (scalar) or 3 dims (vector), Current shape of field is:",field.shape)

    def ComputeConvectiveContribution(self, N, grad_u, a, out):
        """
        Compute (w,a·∇u) although with the definition we employ for ∇u this is actually (w,∇u·a)

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
        if grad_u.ndim == 3: #gradient of a vector function
            xp.einsum('i,ekl,el->eik', N, grad_u, a, out=out, optimize=True)
        elif grad_u.ndim == 2: #gradient of a scalar function
            xp.einsum('i,el,el->ei', N, grad_u, a, out=out, optimize=True)
        else:
            raise ValueError("grad_u must have 2 dims (scalar) or 3 dims (vector)")

    def ComputeMomentumStabilization(self, N, DN, a, u_elemental, Pi_elemental, out, a_DN, PiContrib):
        ##TODO: avoid temporaries!
        xp.einsum("el,eil->ei",a,DN, out=a_DN, optimize=True) #TODO: reuse an auxiliary array
        xp.einsum("eI,eJ,eJk->eIk",a_DN,a_DN,u_elemental,out=out, optimize=True)
        xp.einsum("eI,J,eJk->eIk",a_DN,N,Pi_elemental, out=PiContrib, optimize=True)
        out -= PiContrib

    def ComputeDivDivStabilization(self, DN, Pi_elemental, out):
        ##TODO: avoid temporaries!
        pass

    def ComputePressureStabilization_ProjectionTerm(self, N, DN, Pi_press_el, out):
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
        xp.einsum("eIk,J,eJk->eI",DN,N,Pi_press_el,out=out)

    #REMARK: this function has to be run on the cpu AND NOT ON THE GPU - it will return a gpu matrix if that is the case
    def AllocateScalarMatrix(self,conn : np.ndarray): 
        #TODO: make it efficient
        size = np.max(conn)+1
        Agraph = KM.SparseContiguousRowGraph(size)
        for e in range(conn.shape[0]):
            Agraph.AddEntries(conn[e])
        Agraph.Finalize()

        #assembling matrix
        A = KM.CsrMatrix(Agraph)

        if USE_CUPY:
            assembly_indices = A.GetEquationIdCsrIndices(conn) 
            A_cu = sparse.csr_matrix(
                (xp.asarray(A.value_data(), dtype=PRECISION),
                xp.asarray(A.index2_data()),
                xp.asarray(A.index1_data())),
                shape=(A.Size1(), A.Size2())
            )   
            return A_cu,assembly_indices
        else:
            assembly_indices = A.GetEquationIdCsrIndices(conn) 
            return A,assembly_indices

    def AssembleScalarMatrix(self,LHSel,conn,A : KM.CsrMatrix):
        A.SetValue(0.0)
        A.BeginAssemble()
        for e in range(conn.shape[0]):
            A.Assemble(LHSel[e],conn[e])
        A.FinalizeAssemble()

    def AssembleScalarMatrixByCSRIndices(self,LHSel,csr_indices,A):
        if USE_CUPY:
            A.data.fill(0.0)
            self.AssembleVector(csr_indices.ravel(), LHSel.ravel(), A.data)
        else:
            A.value_data().fill(0.0)
            self.AssembleVector(csr_indices.ravel(), LHSel.ravel(), A.value_data())

    def ApplyHomogeneousDirichlet(self, fixed_values, LHS, RHS, diag_value=1.0):
        """
        Applies homogeneous Dirichlet BCs (value=0) to a CSR matrix and RHS.
        Maintains symmetry by zeroing both rows and columns of fixed indices.
        """
        

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
        # LHS.indptr tells us where rows start; np.diff tells us how many entries per row
        diff_array = xp.diff(indptr).astype(xp.int64)
        row_map = xp.repeat(xp.arange(n), diff_array)
        
        # 3. Create masks for entries in fixed rows and fixed columns
        mask_row = is_fixed[row_map]
        mask_col = is_fixed[indices]
        
        # 4. Zero out the "Symmetry Cross" (all entries in fixed rows OR columns)
        # This is the most expensive step, done in a single vectorized pass
        data[mask_row | mask_col] = 0.0
        
        # 5. Restore the diagonal for the fixed rows
        # A diagonal entry is where row_index == column_index
        mask_diag = (row_map == indices) & mask_row
        data[mask_diag] = diag_value
        
        # 6. Set RHS to 0 at fixed indices
        RHS[fixed_values] = 0.0
        
        return LHS, RHS