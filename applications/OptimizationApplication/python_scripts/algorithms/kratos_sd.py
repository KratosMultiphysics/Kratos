# Kratos Core and Apps
import KratosMultiphysics as KM
import KratosMultiphysics.OptimizationApplication as KOA
from KratosMultiphysics import Parameters, Logger

# Additional imports
from KratosMultiphysics.OptimizationApplication.algorithms.algorithm_base import OptimizationAlgorithm
from KratosMultiphysics.OptimizationApplication.utilities.timer import Timer



# ==============================================================================
class AlgorithmSteepestDescent(OptimizationAlgorithm):
    # --------------------------------------------------------------------------
    def __init__(self):
        self.objective # Standaratize Objective Response 
        self.des_var
        self.MC
        self.grad_f
        self.search_direction
        self.des_var_update
        self.alfa
        self.converged = False


    # --------------------------------------------------------------------------
    def Initialize(self):
        super().InitializeOptimizationLoop()
        self.opt_algorithm.Initialize()
        self._InitializeCSVLogger()
        self.des_var = self.MC.GetInitialDesignVariables()

    # --------------------------------------------------------------------------

    def ComputeSearchDirection(self):
        self.search_direction = - self.gradients

    def QuadApproximation(self, f, test_f, grad_f):
        alfa = 1 +2 
        return alfa

    def LineSearch(self):
        test_des_var_update = self.alfa * self.search_direction
        test_f = self.objective.ComputeResponseValue(self.des_var + test_des_var_update, self.MC)
        return self.QuadApproximation(self.value_f, test_f, self.grad_f)

    def ComputeDesignVariableUpdate(self):
        self.alfa = self.LineSearch()
        self.des_var_update = self.alfa * self.search_direction

    def RunOptimizationLoop(self):

        timer = Timer()
        timer.StartTimer()  

        self.Initialize()

        while not self.converged:

            self.value_f = self.objective.ComputeResponseValue(self.des_var, self.MC)
            self.objective.ComputeResponseGradients(self.des_var, self.grad_f, self.MC)
            self.ComputeSearchDirection()
            self.ComputeDesignVariableUpdate()
            self.des_var += self.des_var_update

            self.CheckConvergence()

        self.Finalize()

    # --------------------------------------------------------------------------
    def Finalize(self):
        pass

# ==============================================================================
