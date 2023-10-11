import numpy as np
from scipy.linalg import qr as scipy_qr

class DiscreteEmpiricalInterpolationMethod():
    def __init__(self, DEIM_tolerance=1e-6):
        self.DEIM_tolerance = DEIM_tolerance
        self.P = None  # Interpolation indices
        self.P_matrix = None  # Matrix constructed from self.P indices
        self.W = None  # Weight matrix

    def SetUp(self, U_k):
        """Set up the Q-DEIM algorithm with the given matrix U_k."""
        self.U_k = U_k

    def Initialize(self):
        """Initialize the Q-DEIM algorithm by selecting the first index."""
        # Choose the first index by finding the row with the maximum norm
        norms = np.linalg.norm(self.U_k, axis=1)
        index = np.argmax(norms)
        self.P = np.array([index])

    def DEIM(self):
        """Calculate the DEIM indices using the standard DEIM algorithm."""
        m, n = self.U_k.shape

        # First DEIM point
        nmax = np.argmax(np.abs(self.U_k[:, 0]))
        XI_m = self.U_k[:, 0].reshape(m, 1)
        z = np.zeros((m, 1))
        self.P_matrix = np.copy(z)
        self.P_matrix[nmax, 0] = 1
        self.P = [nmax]

        # DEIM points 2 to n
        for jj in range(1, n):
            c = np.linalg.solve(self.P_matrix.T @ XI_m, self.P_matrix.T @ self.U_k[:, jj].reshape(m, 1))
            res = self.U_k[:, jj].reshape(m, 1) - XI_m @ c
            nmax = np.argmax(np.abs(res))
            XI_m = np.concatenate((XI_m, self.U_k[:, jj].reshape(m, 1)), axis=1)
            self.P_matrix = np.concatenate((self.P_matrix, z), axis=1)
            self.P_matrix[nmax, jj] = 1
            self.P.append(nmax)
        self.P = np.sort(self.P)  # Ensure self.P is sorted
        self.ConstructInterpolationIndicesMatrix()
        self.ComputeWeightsAndSave()

    def QDEIM(self):
        """Calculate the DEIM indices using the Q-DEIM algorithm."""
        m, n = self.U_k.shape

        # QR factorization with column pivoting
        Q, R, pivot = scipy_qr(self.U_k, pivoting=True)
        
        # Construct the QDEIM matrix P_matrix
        self.P_matrix = np.zeros((m, n))
        for jj in range(n):
            self.P_matrix[pivot[jj], jj] = 1

        # Extract the list of selected indices from the pivot
        self.P = pivot.tolist()
        self.P = np.sort(self.P)  # Ensure self.P is sorted
        self.ConstructInterpolationIndicesMatrix()
        self.ComputeWeightsAndSave()

    def Calculate(self, method="DEIM"):
        """Calculate the DEIM indices using the specified method."""
        if method == "DEIM":
            self.DEIM()
        elif method == "QDEIM":
            self.QDEIM()
        else:
            raise ValueError("Invalid method. Choose either 'DEIM' or 'QDEIM'.")

    def Run(self):
        """Run the Q-DEIM algorithm."""
        # self.Initialize()
        self.Calculate()

    def ComputeWeightsAndSave(self):
        """Compute InterpolationWeights matrix and save InterpolationWeights and InterpolationIndicesMatrix."""
        # Select rows from U_k using self.P
        selected_U_k = self.U_k[self.P, :]
        
        # Compute InterpolationWeights
        self.InterpolationWeights = self.U_k @ np.linalg.pinv(selected_U_k.T)
        
        # Save InterpolationWeights and self.P as numpy arrays
        np.save('InterpolationWeights.npy', self.InterpolationWeights)
        np.save('DEIM_InterpolationIndices.npy', self.P)

    def ConstructInterpolationIndicesMatrix(self):
        """Construct the InterpolationIndicesMatrix from the self.P indices."""
        m, n = self.U_k.shape
        self.InterpolationIndicesMatrix = np.zeros((m, len(self.P)))
        self.InterpolationIndicesMatrix[self.P, np.arange(len(self.P))] = 1

    # def ComputeWeightsAndSave(self):
    #     """Compute InterpolationWeights matrix and save InterpolationWeights and InterpolationIndicesMatrix."""
    #     # Compute InterpolationWeights
    #     self.InterpolationWeights = self.U_k @ np.linalg.pinv(self.InterpolationIndicesMatrix.T @ self.U_k)
        
    #     # Save InterpolationWeights and InterpolationIndicesMatrix as numpy arrays
    #     np.save('InterpolationWeights.npy', self.InterpolationWeights)
    #     np.save('InterpolationIndicesMatrix.npy', self.InterpolationIndicesMatrix)




