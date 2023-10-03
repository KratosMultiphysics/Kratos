import numpy as np

class QDiscreteEmpiricalInterpolationMethod():
    def __init__(self, DEIM_tolerance=1e-6):
        self.DEIM_tolerance = DEIM_tolerance
        self.P = None  # Interpolation indices

    def SetUp(self, U_k):
        self.U_k = U_k

    def Initialize(self):
        # Choose the first index by finding the row with the maximum norm
        norms = np.linalg.norm(self.U_k, axis=1)
        index = np.argmax(norms)
        self.P = np.array([index])

    def Calculate(self):
        m, n = self.U_k.shape
        P_matrix = np.zeros((len(self.P), m))
        P_matrix[np.arange(len(self.P)), self.P] = 1  # Use advanced indexing
        B = np.linalg.pinv(self.U_k[self.P, :])
        for _ in range(n - len(self.P)):
            r = P_matrix @ self.U_k @ B - np.eye(len(self.P))
            q = np.argmax(np.abs(r))
            i, j = divmod(q, len(self.P))
            self.P = np.append(self.P, self.P[i])
            P_matrix = np.zeros((len(self.P), m))
            P_matrix[np.arange(len(self.P)), self.P] = 1  # Use advanced indexing
            B = np.linalg.pinv(self.U_k[self.P, :])

    def QDEIM(self):
        # QR factorization with column pivoting
        Q, R, P = np.linalg.qr(self.U_k, pivoting=True)
        self.P = P[:self.U_k.shape[1]]

    def Run(self):
        self.Initialize()
        self.Calculate()
        self.QDEIM()
        return self.P  # Return the selected indices

# Usage
U_k = np.array([[1, 2], [3, 4], [5, 6]])  # This is an example; replace with your truncated basis
qdeim = QDiscreteEmpiricalInterpolationMethod(DEIM_tolerance=0.01)
qdeim.SetUp(U_k)
selected_indices = qdeim.Run()
print(f"Selected row indices for Q-DEIM: {selected_indices}")




