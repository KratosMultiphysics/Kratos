import numpy as np
import cupy as xp
import scipy.sparse as sp
import cupyx.scipy.sparse as csp
import cupyx.scipy.sparse.linalg as csp_linalg

sparse = csp
sparse_linalg = csp_linalg
from cupyx.scipy.sparse.linalg import LinearOperator, cg
import pyamg
import KratosMultiphysics as KM


class GraphBasedSA_AMG:
    def __init__(self, A_pattern_gpu, max_coarse=300, smoothed_aggregation = False):
        self.smoothed_aggregation = smoothed_aggregation
        if(type(A_pattern_gpu) != KM.CsrMatrix):
            self.shape = A_pattern_gpu.shape
            self.dtype = A_pattern_gpu.dtype
            self.max_coarse = max_coarse

            # ---------------------------------------------------------
            # CPU PHASE: Purely structural graph analysis
            # ---------------------------------------------------------
            indptr_cpu = A_pattern_gpu.indptr.get()
            indices_cpu = A_pattern_gpu.indices.get()
            dummy_data = A_pattern_gpu.data.get()
        else:
            self.shape = xp.array([A_pattern_gpu.size1(), A_pattern_gpu.size2()], dtype=xp.int32)
            self.dtype = A_pattern_gpu.dtype
            self.max_coarse = max_coarse

            indptr_cpu = A_pattern_gpu.index1_data()
            indices_cpu = A_pattern_gpu.index2_data()
            dummy_data = A_pattern_gpu.value_data()

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


    def estimate_lambda_max_gpu(self,A, invD, n_iter=5):
        """
        Estimates lambda_max(D^-1 A) on GPU using Lanczos iteration on 
        the symmetric operator S = D^-1/2 * A * D^-1/2.
        """
        dtype = A.dtype
        n = A.shape[0]
        
        # 1. Symmetric scaling factor D^-1/2
        invD_sqrt = xp.sqrt(invD)
        
        # 2. Seed initial normalized random vector on GPU
        rng = xp.random.default_rng(seed=42)
        v = rng.standard_normal(n, dtype=dtype)
        v /= xp.linalg.norm(v)
        
        v_prev = xp.zeros_like(v)
        beta = 0.0
        
        alphas = []
        betas = []
        
        # 3. Lanczos Iteration Loop
        for _ in range(n_iter):
            # Apply S * v = D^-1/2 * (A * (D^-1/2 * v))
            w = invD_sqrt * (A @ (invD_sqrt * v))
            
            # Compute Ritz coefficient alpha = v^T * S * v
            alpha = float(xp.dot(v, w))
            alphas.append(alpha)
            
            # Orthogonalize against previous vectors
            w -= alpha * v + beta * v_prev
            beta = float(xp.linalg.norm(w))
            
            if beta < 1e-12:
                break
                
            betas.append(beta)
            v_prev = v
            v = w / beta
            
        # 4. Construct small tridiagonal matrix T on CPU (size k x k)
        k = len(alphas)
        T = np.zeros((k, k), dtype=np.float64)
        
        # Fill main diagonal
        np.fill_diagonal(T, alphas)
        
        # Fill off-diagonals safely (a k x k matrix has at most k-1 off-diagonals)
        num_off_diags = min(k - 1, len(betas))
        for i in range(num_off_diags):
            T[i, i + 1] = betas[i]
            T[i + 1, i] = betas[i]
                
        # Compute maximum eigenvalue of the tiny k x k matrix
        eigvals = np.linalg.eigvalsh(T)
        lambda_max = float(np.max(eigvals))
        
        return max(lambda_max, 1.0)

    def update_matrix_values(self, A_new, lambda_estima_safety_factor=1.1, lambda_estimate_iterations=5):
        """
        Phase 2: GPU only - Fast rebuild of P, R, and A_coarse.
        """
        self.levels = []
        A_current = A_new #.tocsr()

        # Lock omegas to prevent Float64 upcasting
        #omega_prolong = xp.float32(4.0 / 3.0) if self.dtype == xp.float32 else 4.0 / 3.0
        omega_prolong = self.dtype.type(4.0/3.0)

        for lvl_idx, T in enumerate(self.T_operators_gpu):
            lvl_dict = {'A': A_current}

            # Let CuPy find and extract the diagonal natively on the GPU
            diag = A_current.diagonal()

            # Protect against near-zeros causing float32 explosions
            #diag[xp.abs(diag) < 1e-7] = 1.0
            invD = xp.reciprocal(diag) #
            #invD = self.dtype(1.0) / diag #TODO: check precision!
            lvl_dict['invD'] = invD

            if self.smoothed_aggregation:
                # Build Smoothed Prolongator using @
                A_T = A_current @ T
                scaled_A_T = sparse.spdiags(invD, 0, A_current.shape[0], A_current.shape[1]) * A_T

                P = T - (scaled_A_T * omega_prolong)
            else: #non smoothed version
                P = T 
            lvl_dict['P'] = P
            R = P.T

            # ------------------------------------------------------------------
            # Compute LEVEL-SPECIFIC optimal Jacobi weight on GPU
            # ------------------------------------------------------------------
            # 5 steps of Power Iteration / Lanczos on current level's (invD * A_current)
            lambda_max_lvl = self.estimate_lambda_max_gpu(A_current, invD, n_iter=lambda_estimate_iterations)
            
            # Apply safety margin (1.05) and compute optimal damping: (4/3) / lambda_max
            safe_lambda = max(lambda_max_lvl, 1.0) * lambda_estima_safety_factor
            lvl_dict['omega_jacobi'] = omega_prolong / safe_lambda
            #print("omega level = ", lvl_dict['omega_jacobi'])

            self.levels.append(lvl_dict)

            # Generate the new structure natively on GPU using @
            A_current = R @ A_current @ P

        # --- Base (Coarsest) Level Handling ---
        base_dict = {'A': A_current}

        if A_current.shape[0] <= self.max_coarse:
            base_dict['is_dense_solve'] = True
            
            # A_dense_64 = A_current.todense().astype(xp.float64)
            # # Use pseudo-inverse to survive singular coarse grids
            # A_inv_64 = xp.linalg.pinv(A_dense_64)
            # base_dict['A_dense_inv'] = A_inv_64.astype(self.dtype)

            A_dense = A_current.todense()
            # Use pseudo-inverse to survive singular coarse grids
            #Safe threshold: machine epsilon * max dimension * 100
            A_inv = xp.linalg.pinv(A_dense, rcond=A_dense.shape[0] * xp.finfo(self.dtype).eps * 10)
            #enforce symmetry
            A_inv = (A_inv + A_inv.T) / 2.
            base_dict['A_dense_inv'] = A_inv

        else:
            base_dict['is_dense_solve'] = False

            # Native extraction for the base level
            diag = A_current.diagonal()
            diag[xp.abs(diag) < 1e-7] = 1.0
            base_dict['invD'] = self.dtype.type(1.0) / diag

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
                omega_jacobi = lvl["omega_jacobi"] #  self.dtype(0.67)
                for _ in range(20):
                    r = b - (A @ x)
                    r *= lvl['invD']
                    r *= omega_jacobi
                    x += r
                return x

        P = lvl['P']
        invD = lvl['invD']
        omega = lvl["omega_jacobi"] # self.dtype.type(0.6) #self.dtype.type(0.67)
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

