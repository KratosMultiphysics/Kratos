import numpy as np
from scipy.sparse.linalg import gmres, LinearOperator, aslinearoperator


from KratosMultiphysics.CoSimulationApplication.coupled_solvers.gauss_seidel import CoupledSolverGaussSeidel
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools
cs_data_structure = cs_tools.cs_data_structure


def Create(parameters):
    return CoupledSolverIBQN(parameters)

class CoupledSolverIBQN(CoupledSolverGaussSeidel):
    def __init__(self, parameters):
        """This coupled solver is not working properly yet!"""

        super().__init__(parameters)

        self.model_f = cs_tools.CreateInstance(self.parameters["settings"]["model"])
        self.model_s = cs_tools.CreateInstance(self.parameters["settings"]["model"])
        self.omega = self.settings["omega"].GetDouble()
        self.rtol = self.settings["relative_tolerance_gmres"].GetDouble()

        self.y = None
        self.xtemp = None
        self.ytemp = None
        self.u = None
        self.w = None

    def Initialize(self):
        super().Initialize()

        self.y = self.solver_wrappers[0].GetInterfaceOutput()
        self.xtemp = self.x.deepcopy()
        self.ytemp = self.y.deepcopy()
        self.u = self.x.GetNumpyArray().shape[0]
        self.w = self.y.GetNumpyArray().shape[0]
        self.model_f.size = self.u
        self.model_s.size = self.w
        models = [self.model_f, self.model_s]
        for model in models:
            model.Initialize()
            self.components += [model]

    def lop_f(self, x):
        xtemp = self.xtemp.deepcopy()
        xtemp.SetNumpyArray(x.flatten())
        return self.model_f.Predict(xtemp).GetNumpyArray()

    def lop_s(self, y):
        ytemp = self.ytemp.deepcopy()
        ytemp.SetNumpyArray(y.flatten())
        return self.model_s.Predict(ytemp).GetNumpyArray()

    def resx(self, rk):
        print("x " + str(rk))

    def resy(self, rk):
        print("y " + str(rk))

    def SolveSolutionStep(self):
        iu = aslinearoperator(np.identity(self.u))
        iw = aslinearoperator(np.identity(self.w))
        dx = self.x.deepcopy()
        dy = self.x.deepcopy()
        # Initial values
        self.x = self.predictor.Predict(self.x)
        yt = self.solver_wrappers[0].SolveSolutionStep(self.x)
        print('x = ' + str(np.linalg.norm(self.x.GetNumpyArray())), 'yt = ' + str(np.linalg.norm(yt.GetNumpyArray())))
        self.model_f.Add(self.x, yt)
        self.y = yt
        xt = self.solver_wrappers[1].SolveSolutionStep(self.y)
        print('x = ' + str(np.linalg.norm(self.x.GetNumpyArray())), 'yt = ' + str(np.linalg.norm(yt.GetNumpyArray())))
        print('y = ' + str(np.linalg.norm(self.y.GetNumpyArray())), 'xt = ' + str(np.linalg.norm(xt.GetNumpyArray())))
        self.model_s.Add(self.y, xt)
        r = xt - self.x
        self.convergence_criterion.Update(r)
        # Coupling iteration loop
        while not self.convergence_criterion.IsSatisfied():
            if not self.model_s.IsReady() or not self.model_f.IsReady:
                dx = self.omega * r
                print("not ready 1")
            else:
                mf = LinearOperator((self.w, self.u), self.lop_f)
                ms = LinearOperator((self.u, self.w), self.lop_s)
                a = iu - ms @ mf
                b = (xt - self.x).GetNumpyArray() + ms @ (yt - self.y).GetNumpyArray()
                dx_sol, exitcode = gmres(a, b, tol=self.rtol)
                print("x max = " + str(np.max(a * dx_sol - b)), "norm = " + str(np.linalg.norm(a * dx_sol - b)))
                if exitcode != 0:
                    raise RuntimeError("GMRES failed")
                dx.SetNumpyArray(dx_sol)
            self.x += dx
            yt = self.solver_wrappers[0].SolveSolutionStep(self.x)
            print('y = ' + str(np.linalg.norm(self.y.GetNumpyArray())), 'xt = ' + str(np.linalg.norm(xt.GetNumpyArray())))
            print('x = ' + str(np.linalg.norm(self.x.GetNumpyArray())), 'yt = ' + str(np.linalg.norm(yt.GetNumpyArray())))
            self.model_f.Add(self.x, yt)
            if not self.model_s.IsReady() or not self.model_f.IsReady:
                dy = yt - self.y
                print("not ready 2")
            else:
                a = iw - mf @ ms
                b = (yt - self.y).GetNumpyArray() + mf @ (xt - self.x).GetNumpyArray()
                dy_sol, exitcode = gmres(a, b, tol=self.rtol)
                print("y max = " + str(np.max(a * dy_sol - b)), "norm = " + str(np.linalg.norm(a * dy_sol - b)))
                if exitcode != 0:
                    raise RuntimeError("GMRES failed")
                dy.SetNumpyArray(dy_sol)
            self.y += dy
            xt = self.solver_wrappers[1].SolveSolutionStep(self.y)
            print('x = ' + str(np.linalg.norm(self.x.GetNumpyArray())), 'yt = ' + str(np.linalg.norm(yt.GetNumpyArray())))
            print('y = ' + str(np.linalg.norm(self.y.GetNumpyArray())), 'xt = ' + str(np.linalg.norm(xt.GetNumpyArray())))
            self.model_s.Add(self.y, xt)
            r = xt - self.x
            self.convergence_criterion.Update(r)


