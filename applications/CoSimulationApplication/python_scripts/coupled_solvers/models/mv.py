import numpy as np

from KratosMultiphysics.CoSimulationApplication.co_simulation_component import CoSimulationComponent
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools
cs_data_structure = cs_tools.cs_data_structure


def Create(parameters):
    return ModelMV(parameters)


class ModelMV(CoSimulationComponent):
    def __init__(self, parameters):
        super().__init__()

        self.settings = parameters["settings"]
        self.min_significant = self.settings["min_significant"].GetDouble()

        self.size = None
        self.added = False
        self.rref = None
        self.xtref = None
        self.v = None
        self.w = None
        self.ncurr = None
        self.nprev = None

    def Initialize(self):
        super().Initialize()

        self.v = np.empty((self.size, 0))
        self.w = np.empty((self.size, 0))
        self.nprev = np.zeros((self.size, self.size))

    def Filter(self):
        if self.v.shape[1] == 0:
            raise RuntimeError("No information to filter")
        # Remove columns resulting in small diagonal elements in R
        singular = True
        while singular and self.v.shape[1]:
            rr = np.linalg.qr(self.v, mode='r')
            diag = np.diagonal(rr)
            m = min(abs(diag))
            if m < self.min_significant:
                i = np.argmin(abs(diag))
                print("Removing column " + str(i) + ": " + str(m) + " < minsignificant")
                self.v = np.delete(self.v, i, 1)
                self.w = np.delete(self.w, i, 1)
            else:
                singular = False
        # Remove columns if number of columns exceeds number of rows
        if self.v.shape[0] < self.v.shape[1]:
            self.v = np.delete(self.v, -1, 1)
            self.w = np.delete(self.w, -1, 1)

    def Predict(self, r_in):
        r = r_in.GetNumpyArray().reshape(-1, 1)
        if self.ncurr is None:
            raise RuntimeError("No information to predict")
        # Approximation for the inverse of the Jacobian from a multiple vector model
        dr = -r
        dxt = self.ncurr @ dr
        dxt_out = r_in.deepcopy()
        dxt_out.SetNumpyArray(dxt.flatten())
        return dxt_out

    def Add(self, r_in, xt_in):
        r = r_in.GetNumpyArray().reshape(-1, 1)
        xt = xt_in.GetNumpyArray().reshape(-1, 1)
        if self.added:
            dr = r - self.rref
            dxt = xt - self.xtref
            # Update V and W matrices
            self.v = np.hstack((dr, self.v))
            self.w = np.hstack((dxt, self.w))
            self.Filter()
            # Update of the matrix N
            self.ncurr = self.nprev + (self.w - self.nprev @ self.v) @ np.linalg.inv(self.v.T @ self.v) @ self.v.T
        else:
            self.added = True
        self.rref = r
        self.xtref = xt

    def IsReady(self):
        return self.ncurr is not None

    def InitializeSolutionStep(self):
        super().InitializeSolutionStep()

        self.rref = None
        self.xtref = None
        self.v = np.empty((self.size, 0))
        self.w = np.empty((self.size, 0))
        self.added = False

    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()

        self.nprev = self.ncurr
