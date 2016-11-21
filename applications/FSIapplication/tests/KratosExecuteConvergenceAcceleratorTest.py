from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

from KratosMultiphysics import *
from KratosMultiphysics.FSIApplication import *

import process_factory
import KratosMultiphysics.KratosUnittest as KratosUnittest

class KratosExecuteConvergenceAcceleratorTest(KratosUnittest.TestCase):

    def __init__(self, ProjectParameters):
        
        self.vector_space = UblasSparseSpace()
        
        self.settings = Parameters("{}")
        self.settings.AddValue("coupling_solver_settings",ProjectParameters["coupling_solver_settings"]["solver_settings"])
        coupling_utility_parameters = self.settings["coupling_solver_settings"]["coupling_strategy"]
        
        # Construct the coupling partitioned strategy
        import convergence_accelerator_factory     
        self.coupling_utility = convergence_accelerator_factory.CreateConvergenceAccelerator(coupling_utility_parameters)
        print("* Coupling strategy constructed.")
        
    def Solve(self):
        
        tol = self.settings["coupling_solver_settings"]["nl_tol"].GetDouble()
        max_it = self.settings["coupling_solver_settings"]["nl_max_it"].GetInt()
        recursive_test = self.settings["coupling_solver_settings"]["recursive_test"].GetBool()
        
        self.coupling_utility.Initialize()
        
        if recursive_test == True:
            recursive_steps = 10
        else:
            recursive_steps = 1
        
        step = 1
        res_norm = 1.0

        while (step <= recursive_steps):
            
            nl_it = 1
            convergence = False
            
            # x_guess initialization
            x_guess = Vector(6)
            for i in range(0,len(x_guess)):
                x_guess[i] = step
            
            self.coupling_utility.InitializeSolutionStep()
            
            while (nl_it <= max_it):
                
                residual = ComputeResidual(x_guess)
                res_norm = self.vector_space.TwoNorm(residual)
                
                print("Iteration: ", nl_it," residual norm: ", res_norm)
                            
                if res_norm > tol:
                    self.coupling_utility.UpdateSolution(residual, x_guess)
                    self.coupling_utility.FinalizeNonLinearIteration()
                    nl_it += 1
                    
                else:
                    self.coupling_utility.FinalizeSolutionStep()
                    convergence = True
                    break
                
            step += 1
        
        # Check the obtained solution
        expected_RHS = Vector(6)
        expected_RHS[0] = -3.49458887
        expected_RHS[1] = 11.98381557
        expected_RHS[2] = -9.42225393
        expected_RHS[3] = 3.79688473 
        expected_RHS[4] = 0.72338935
        expected_RHS[5] = 4.36818068
        
        obtained_RHS = ComputeExpectedRHS(x_guess)
        
        if convergence == True:
            for i in range(0,len(expected_RHS)):
                self.assertAlmostEqual(obtained_RHS[i],expected_RHS[i])        
        
            
def ComputeResidual(x_guess):
    
    res = Vector(6)
    
    res[0] = (-3.49458887) - (0.98071655*x_guess[0] + 0.19229810*x_guess[1] + 0.30490177*x_guess[2] + 0.23445587*x_guess[3] + 0.01612071*x_guess[4] + 0.70824327*x_guess[5])
    res[1] = (11.98381557) - (0.17029604*x_guess[0] + 0.15202212*x_guess[1] + 0.49687295*x_guess[2] + 0.83274282*x_guess[3] + 0.32594298*x_guess[4] + 0.47098796*x_guess[5])
    res[2] = (-9.42225393) - (0.87776675*x_guess[0] + 0.55563239*x_guess[1] + 0.48020994*x_guess[2] + 0.01211062*x_guess[3] + 0.77612792*x_guess[4] + 0.99036026*x_guess[5])
    res[3] = (3.796884730) - (0.88236065*x_guess[0] + 0.87559386*x_guess[1] + 0.74122791*x_guess[2] + 0.91010412*x_guess[3] + 0.25651983*x_guess[4] + 0.44262038*x_guess[5])
    res[4] = (0.723389350) - (0.53553728*x_guess[0] + 0.44427695*x_guess[1] + 0.22965224*x_guess[2] + 0.51367362*x_guess[3] + 0.95841073*x_guess[4] + 0.93117203*x_guess[5])
    res[5] = (4.368180680) - (0.74119728*x_guess[0] + 0.80883251*x_guess[1] + 0.06562307*x_guess[2] + 0.14071936*x_guess[3] + 0.28756793*x_guess[4] + 0.63488182*x_guess[5])
    
    return res
    
def ComputeExpectedRHS(x_guess):
    
    RHS = Vector(6)

    RHS[0] = (0.98071655*x_guess[0] + 0.19229810*x_guess[1] + 0.30490177*x_guess[2] + 0.23445587*x_guess[3] + 0.01612071*x_guess[4] + 0.70824327*x_guess[5])
    RHS[1] = (0.17029604*x_guess[0] + 0.15202212*x_guess[1] + 0.49687295*x_guess[2] + 0.83274282*x_guess[3] + 0.32594298*x_guess[4] + 0.47098796*x_guess[5])
    RHS[2] = (0.87776675*x_guess[0] + 0.55563239*x_guess[1] + 0.48020994*x_guess[2] + 0.01211062*x_guess[3] + 0.77612792*x_guess[4] + 0.99036026*x_guess[5])
    RHS[3] = (0.88236065*x_guess[0] + 0.87559386*x_guess[1] + 0.74122791*x_guess[2] + 0.91010412*x_guess[3] + 0.25651983*x_guess[4] + 0.44262038*x_guess[5])
    RHS[4] = (0.53553728*x_guess[0] + 0.44427695*x_guess[1] + 0.22965224*x_guess[2] + 0.51367362*x_guess[3] + 0.95841073*x_guess[4] + 0.93117203*x_guess[5])
    RHS[5] = (0.74119728*x_guess[0] + 0.80883251*x_guess[1] + 0.06562307*x_guess[2] + 0.14071936*x_guess[3] + 0.28756793*x_guess[4] + 0.63488182*x_guess[5])

    return RHS
