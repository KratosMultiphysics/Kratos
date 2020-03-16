import numpy as np
from scipy.sparse.linalg import gmres, LinearOperator, aslinearoperator


from KratosMultiphysics.CoSimulationApplication.coupled_solvers.gauss_seidel import CoupledSolverGaussSeidel
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools
cs_data_structure = cs_tools.cs_data_structure


def Create(parameters):
    return CoupledSolverIBQN(parameters)

class CoupledSolverIBQN(CoupledSolverGaussSeidel):
    def __init__(self, parameters):
        super().__init__(parameters)

        self.model_f = cs_tools.CreateInstance(self.parameters["settings"]["model_f"])
        self.model_s = cs_tools.CreateInstance(self.parameters["settings"]["model_s"])
        self.omega = self.settings["omega"].GetDouble()
        self.atol = self.settings["absolute_tolerance_gmres"].GetDouble()
        self.rtol = self.settings["relative_tolerance_gmres"].GetDouble()

        self.y = None
        self.xtemp = None
        self.ytemp = None
        self.u = None
        self.w = None
        self.ready = None

    def Initialize(self):
        super().Initialize()

        self.y = self.solver_wrappers[0].GetInterfaceOutput()
        self.xtemp = self.x.deepcopy()
        self.ytemp = self.y.deepcopy()
        self.u = self.x.GetNumpyArray().shape[0]
        self.w = self.y.GetNumpyArray().shape[0]
        self.ready = False
        self.model_f.size_in = self.model_s.size_out = self.u
        self.model_f.size_out = self.model_s.size_in = self.w
        self.model_f.out = self.y.deepcopy()
        self.model_s.out = self.x.deepcopy()
        models = [self.model_f, self.model_s]
        for model in models:
            model.Initialize()
            self.components += [model]

    def lop_f(self, x):
        self.xtemp.SetNumpyArray(x.flatten())
        return self.model_f.Predict(self.xtemp).GetNumpyArray()

    def lop_s(self, y):
        self.ytemp.SetNumpyArray(y.flatten())
        return self.model_s.Predict(self.ytemp).GetNumpyArray()

    def SolveSolutionStep(self):
        iu = aslinearoperator(np.identity(self.u))
        iw = aslinearoperator(np.identity(self.w))
        dx = self.x.deepcopy()
        dy = self.y.deepcopy()
        # Initial value
        self.x = self.predictor.Predict(self.x)
        # First coupling iteration
        yt = self.solver_wrappers[0].SolveSolutionStep(self.x)
        self.model_f.Add(self.x, yt)
        self.y = yt
        xt = self.solver_wrappers[1].SolveSolutionStep(self.y)
        self.model_s.Add(self.y, xt)
        r = xt - self.x
        self.FinalizeIteration(r)
        # Coupling iteration loop
        while not self.convergence_criterion.IsSatisfied():
            if not self.model_s.IsReady() or not self.model_f.IsReady:
                dx = self.omega * r
            else:
                mf = LinearOperator((self.w, self.u), self.lop_f)
                ms = LinearOperator((self.u, self.w), self.lop_s)
                a = iu - ms @ mf
                b = (xt - self.x).GetNumpyArray() + ms @ (yt - self.y).GetNumpyArray()
                dx_sol, exitcode = gmres(a, b, tol=self.rtol, atol=self.atol)
                if exitcode != 0:
                    raise RuntimeError("GMRES failed")
                dx.SetNumpyArray(dx_sol)
            self.x += dx
            yt = self.solver_wrappers[0].SolveSolutionStep(self.x)
            self.model_f.Add(self.x, yt)
            if not self.model_s.IsReady() or not self.model_f.IsReady:
                dy = yt - self.y
            else:
                a = iw - mf @ ms
                b = (yt - self.y).GetNumpyArray() + mf @ (xt - self.x).GetNumpyArray()
                dy_sol, exitcode = gmres(a, b, tol=self.rtol, atol=self.atol)
                if exitcode != 0:
                    raise RuntimeError("GMRES failed")
                dy.SetNumpyArray(dy_sol)
            self.y += dy
            xt = self.solver_wrappers[1].SolveSolutionStep(self.y)
            self.model_s.Add(self.y, xt)
            r = xt - self.x
            self.FinalizeIteration(r)


