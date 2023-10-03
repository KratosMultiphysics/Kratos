import numpy as np

class DiscreteEmpiricalInterpolationMethod():
    def __init__(self, DEIM_tolerance=1e-6, rank=None):
        self.DEIM_tolerance = DEIM_tolerance
        self.rank = rank  # Rank for the reduced basis
        self.P = None  # Interpolation indices

    def SetUp(self, ResidualMatrix):
        self.ResidualMatrix = ResidualMatrix
        # Perform SVD
        U, s, Vt = np.linalg.svd(self.ResidualMatrix)
        # Select the first k columns of U to form a reduced basis U_k
        k = self.rank if self.rank else len(U)
        self.U_k = U[:, :k]

    def Initialize(self):
        # Choose the first index by finding the row with the maximum norm
        norms = np.linalg.norm(self.U_k, axis=1)
        index = np.argmax(norms)
        self.P = np.array([index])

    def Calculate(self):
        m, n = self.ResidualMatrix.shape
        P_matrix = np.zeros((len(self.P), m))
        for i, index in enumerate(self.P):
            P_matrix[i, index] = 1
        B = np.linalg.pinv(self.U_k[self.P, :])
        for _ in range(n - len(self.P)):
            r = P_matrix @ self.U_k @ B - np.eye(len(self.P))
            q = np.argmax(np.abs(r))
            i, j = divmod(q, len(self.P))
            self.P = np.append(self.P, self.P[i])
            P_matrix = np.zeros((len(self.P), m))
            for i, index in enumerate(self.P):
                P_matrix[i, index] = 1
            B = np.linalg.pinv(self.U_k[self.P, :])

    def Run(self):
        self.Initialize()
        self.Calculate()
        return self.P  # Return the selected indices

# Usage
ResidualMatrix = np.array([[1, 2], [3, 4], [5, 6]])
deim = DiscreteEmpiricalInterpolationMethod(DEIM_tolerance=0.01, rank=2)
deim.SetUp(ResidualMatrix)
selected_indices = deim.Run()
print(f"Selected row indices for DEIM: {selected_indices}")



