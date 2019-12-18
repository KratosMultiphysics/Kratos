import numpy as np
from scipy.linalg import solve_triangular

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
        self.q = self.settings["q"].GetDouble()

        self.size = None
        self.added = False
        self.rref = None
        self.xtref = None
        self.v = None
        self.w = None
        self.wprev = []
        self.rrprev = []
        self.qqprev = []

    def Initialize(self):
        super().Initialize()

        self.v = np.empty((self.size, 0))
        self.w = np.empty((self.size, 0))

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
        if not self.v.shape[1] and len(self.wprev) == 0:
            raise RuntimeError("No information to predict")
        # Approximation for the inverse of the Jacobian from a multiple vector model
        dr = -r
        if self.v.shape[1]:
            qq, rr = np.linalg.qr(self.v, mode='reduced')
            b = qq.T @ dr
            c = solve_triangular(rr, b)
            dxt = self.w @ c
        else:
            dxt = np.zeros((self.size, 1))
            qq = np.zeros((self.size, 1))
        i = 0
        while np.linalg.norm(dr) > self.min_significant and i < len(self.wprev):
            dr = dr - qq @ (qq.T @ dr)
            print(str(i) + '   dr = ' + str(np.linalg.norm(dr)))
            b = self.qqprev[i].T @ dr
            c = solve_triangular(self.rrprev[i], b)
            dxt += self.wprev[i] @ c
            qq = self.qqprev[i]
            i += 1
        # Remove insignificant information
        while i < len(self.wprev):
            self.wprev.pop(i)
            self.rrprev.pop(i)
            self.qqprev.pop(i)
            i += 1
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
        else:
            self.added = True
        self.rref = r
        self.xtref = xt

    def IsReady(self):
        return self.v.shape[1] or len(self.wprev) != 0

    def InitializeSolutionStep(self):
        super().InitializeSolutionStep()

        self.rref = None
        self.xtref = None
        self.v = np.empty((self.size, 0))
        self.w = np.empty((self.size, 0))
        self.added = False

    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()

        self.wprev = [self.w] + self.wprev
        qq, rr = np.linalg.qr(self.v, mode='reduced')
        self.rrprev = [rr] + self.rrprev
        self.qqprev = [qq] + self.qqprev
        # Limit number of timesteps reused to q
        if len(self.wprev) > self.q:
            self.wprev.pop()
            self.rrprev.pop()
            self.qqprev.pop()
