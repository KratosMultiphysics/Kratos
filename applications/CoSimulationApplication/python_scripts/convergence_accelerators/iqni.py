import numpy as np
from scipy.linalg import qr, solve_triangular

from KratosMultiphysics.CoSimulationApplication.co_simulation_component import CoSimulationComponent
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools
cs_data_structure = cs_tools.cs_data_structure


def Create(parameters):
    return ConvergenceAcceleratorIQNI(parameters)


class ConvergenceAcceleratorIQNI(CoSimulationComponent):
    def __init__(self, parameters):
        super().__init__()

        self.settings = parameters["settings"]
        self.min_significant = self.settings["min_significant"].GetDouble()
        self.omega = self.settings["omega"].GetDouble()

        self.added = False
        self.rref = np.array([])
        self.xtref = np.array([])
        self.v = np.matrix([])
        self.w = np.matrix([])

    def Predict(self, r_in):
        r = r_in.GetNumpyArray()
        # Remove columns resulting in small diagonal elements in R
        singular = True
        while singular and self.v.shape[1]:
            rr = qr(self.v, mode='r')
            diag = np.diagonal(rr)
            m = min(diag)
            if abs(m) < self.min_significant:
                i = np.argmin(diag)
                self.v = np.delete(self.v, i, 1)
                self.w = np.delete(self.w, i, 1)
                print("Removing columns " + str(i) + ": " + str(abs(m)) + " < minsignificant")
            else:
                singular = False
        # Calculate return value if sufficient data available
        if self.v.shape[1]:
            # Interface Quasi-Newton with approximation for the inverse of the Jacobian from a least-squares model
            qq, rr, *_ = qr(self.v, mode='economic')
            dr = np.reshape(-r, (-1, 1))
            b = qq.T @ dr
            c = solve_triangular(rr, b)
            dx = self.w @ c - dr
            dx = np.asarray(dx).reshape(-1)
        else:
            if self.added:
                dx = self.omega * r
            else:
                raise RuntimeError("No information to predict")
        dx_out = r_in.deepcopy()
        dx_out.SetNumpyArray(dx)
        return dx_out

    def Update(self, x_in, xt_in):
        x = x_in.GetNumpyArray()
        xt = xt_in.GetNumpyArray()
        r = xt - x
        if self.added:
            dr = r - self.rref
            dxt = xt - self.xtref
            dr = np.reshape(dr, (-1, 1))
            dxt = np.reshape(dxt, (-1, 1))
            if self.v.shape[1]:
                self.v = np.concatenate((dr, self.v), 1)
                self.w = np.concatenate((dxt, self.w), 1)
            else:
                self.v = np.matrix(dr)
                self.w = np.matrix(dxt)
        self.rref = r
        self.xtref = np.array(xt)
        self.added = True

    def IsReady(self):
        return self.added

    def InitializeSolutionStep(self):
        super().InitializeSolutionStep()

        self.rref = np.array([])
        self.xtref = np.array([])
        self.v = np.matrix([])
        self.w = np.matrix([])
        self.added = False
