
import KratosMultiphysics
import KratosMultiphysics.FSIApplication as KratosFSI
import KratosMultiphysics.KratosUnittest as KratosUnittest

import math
from KratosMultiphysics.FSIApplication import convergence_accelerator_factory

class ConvergenceAcceleratorTest(KratosUnittest.TestCase):

    # Residual functions
    def f(self,x):
        return 1 - x[0]**2 - x[1]**2

    def g(self,x):
        return x[1] - x[0]

    def ComputeResidual(self,x):
        r = KratosMultiphysics.Vector(2)
        r[0] = self.f(x)
        r[1] = self.g(x)
        return r

    # Specific functions for the Newton-Rahpson acceleration test
    #-Dres/Dx (Note that the minus sign is already included)
    def Jacobian(self,x):
        K = KratosMultiphysics.Matrix(2,2)
        K[0,0] = 2*x[0]
        K[0,1] = 2*x[1]
        K[1,0] = 1.0
        K[1,1] = -1.0
        return K

    def InvertMatrix(self,A):
        Ainv = KratosMultiphysics.Matrix(2,2);
        det = A[0,0]*A[1,1]-A[1,0]*A[0,1]

        Ainv[0,0] = A[1,1]/det
        Ainv[0,1] = -A[0,1]/det
        Ainv[1,0] = -A[1,0]/det
        Ainv[1,1] = A[0,0]/det
        return Ainv

    def prod(self,A,x):
        y = KratosMultiphysics.Vector(2)
        y[0] = A[0,0]*x[0] + A[0,1]*x[1]
        y[1] = A[1,0]*x[0] + A[1,1]*x[1]
        return y

    def norm(self, v):
        return math.sqrt(sum([i**2 for i in v]))

    # Specific residual functions for the MVQN accelerator test
    # Recall that the MVQN accelerator requires a "large enough" problem to avoid becoming ill-conditioned
    def ComputeMVQNResidual(self, x_guess, multiplier = 1.0):

        res = KratosMultiphysics.Vector(6)

        res[0] = multiplier * (-3.49458887) - (0.98071655*x_guess[0] + 0.19229810*x_guess[1] + 0.30490177*x_guess[2] + 0.23445587*x_guess[3] + 0.01612071*x_guess[4] + 0.70824327*x_guess[5])
        res[1] = multiplier * (11.98381557) - (0.17029604*x_guess[0] + 0.15202212*x_guess[1] + 0.49687295*x_guess[2] + 0.83274282*x_guess[3] + 0.32594298*x_guess[4] + 0.47098796*x_guess[5])
        res[2] = multiplier * (-9.42225393) - (0.87776675*x_guess[0] + 0.55563239*x_guess[1] + 0.48020994*x_guess[2] + 0.01211062*x_guess[3] + 0.77612792*x_guess[4] + 0.99036026*x_guess[5])
        res[3] = multiplier * (3.796884730) - (0.88236065*x_guess[0] + 0.87559386*x_guess[1] + 0.74122791*x_guess[2] + 0.91010412*x_guess[3] + 0.25651983*x_guess[4] + 0.44262038*x_guess[5])
        res[4] = multiplier * (0.723389350) - (0.53553728*x_guess[0] + 0.44427695*x_guess[1] + 0.22965224*x_guess[2] + 0.51367362*x_guess[3] + 0.95841073*x_guess[4] + 0.93117203*x_guess[5])
        res[5] = multiplier * (4.368180680) - (0.74119728*x_guess[0] + 0.80883251*x_guess[1] + 0.06562307*x_guess[2] + 0.14071936*x_guess[3] + 0.28756793*x_guess[4] + 0.63488182*x_guess[5])

        return res

    def ComputeMVQNObtainededRHS(self, x_guess):

        RHS = KratosMultiphysics.Vector(6)

        RHS[0] = (0.98071655*x_guess[0] + 0.19229810*x_guess[1] + 0.30490177*x_guess[2] + 0.23445587*x_guess[3] + 0.01612071*x_guess[4] + 0.70824327*x_guess[5])
        RHS[1] = (0.17029604*x_guess[0] + 0.15202212*x_guess[1] + 0.49687295*x_guess[2] + 0.83274282*x_guess[3] + 0.32594298*x_guess[4] + 0.47098796*x_guess[5])
        RHS[2] = (0.87776675*x_guess[0] + 0.55563239*x_guess[1] + 0.48020994*x_guess[2] + 0.01211062*x_guess[3] + 0.77612792*x_guess[4] + 0.99036026*x_guess[5])
        RHS[3] = (0.88236065*x_guess[0] + 0.87559386*x_guess[1] + 0.74122791*x_guess[2] + 0.91010412*x_guess[3] + 0.25651983*x_guess[4] + 0.44262038*x_guess[5])
        RHS[4] = (0.53553728*x_guess[0] + 0.44427695*x_guess[1] + 0.22965224*x_guess[2] + 0.51367362*x_guess[3] + 0.95841073*x_guess[4] + 0.93117203*x_guess[5])
        RHS[5] = (0.74119728*x_guess[0] + 0.80883251*x_guess[1] + 0.06562307*x_guess[2] + 0.14071936*x_guess[3] + 0.28756793*x_guess[4] + 0.63488182*x_guess[5])

        return RHS

    def ComputeMVQNRandomizedSVDResidual(self, x_guess, multiplier):
        A = [
            [3, 2, 0, 0, 2, 3, 2, 5, 0, 0],
            [1,-1, 0, 4, 2,-1, 2, 4, 0, 0],
            [0, 5, 1, 9, 9, 1, 2, 3, 0, 0],
            [1, 3, 1, 3, 5, 7, 3, 2, 0, 0],
            [5, 2,-1, 9, 9, 5, 3, 1, 0, 0],
            [1,-3, 2, 4, 8, 2, 3, 0, 0, 0],
            [5, 2,-1, 9, 9, 5, 4, 3, 0, 0],
            [0, 0,-2,-4, 8, 2, 7, 8, 0, 2],
            [1, 0,-1, 0, 1, 0,-1, 0, 1, 0],
            [0, 0,-2,-4, 8, 2, 7, 8, 2, 3]
        ]
        b = [2, 4,-1, 5, 6,-2, 8,-2, 0, 1]

        res = []
        res = [multiplier*b[i] for i in range(len(b))]
        for i in range(10):
            for j in range(10):
                res[i] -= A[i][j]*x_guess[j]

        return res

    # Aitken accelerator test
    def test_aitken_accelerator(self):

        aitken_settings = KratosMultiphysics.Parameters("""{
                                                            "solver_type"        : "Relaxation",
                                                            "acceleration_type"  : "Aitken",
                                                            "w_0"                : 0.825
                                                           }""")

        print("")
        print("Testing accelerator: ",aitken_settings["solver_type"].GetString())

        # Construct the accelerator strategy
        self.coupling_utility = convergence_accelerator_factory.CreateConvergenceAccelerator(aitken_settings)

        x_guess = KratosMultiphysics.Vector(2)
        x_guess[0] = 0.5
        x_guess[1] = 1.0

        tol = 1e-14
        max_it = 20

        self.coupling_utility.Initialize()

        nl_it = 1
        res_norm = 1.0
        convergence = False

        self.coupling_utility.InitializeSolutionStep()

        while (nl_it <= max_it):

            residual = self.ComputeResidual(x_guess)
            res_norm = math.sqrt(residual[0]**2 + residual[1]**2)

            print("Iteration: ", nl_it," residual norm: ", res_norm)

            if res_norm > tol:
                self.coupling_utility.InitializeNonLinearIteration()
                self.coupling_utility.UpdateSolution(residual, x_guess)
                self.coupling_utility.FinalizeNonLinearIteration()
                nl_it += 1
            else:
                self.coupling_utility.FinalizeSolutionStep()
                convergence = True
                break

        x_final = x_guess

        # Check the obtained solution
        expected_x = KratosMultiphysics.Vector(2)
        expected_x[0] = math.sqrt(2)/2.0
        expected_x[1] = math.sqrt(2)/2.0

        if convergence == True:
            for i in range(0,len(expected_x)):
                self.assertAlmostEqual(expected_x[i],x_final[i])

    # MVQN accelerator test
    def test_mvqn_accelerator(self):

        mvqn_settings = KratosMultiphysics.Parameters("""{
                                                          "solver_type" : "MVQN",
                                                          "w_0"         : 0.825
                                                         }""")

        print("")
        print("Testing accelerator: ",mvqn_settings["solver_type"].GetString())

        # Construct the accelerator strategy
        self.coupling_utility = convergence_accelerator_factory.CreateConvergenceAccelerator(mvqn_settings)

        x_guess = KratosMultiphysics.Vector(6)
        x_guess[0] = -20.0
        x_guess[1] =  20.0
        x_guess[2] =  -5.0
        x_guess[3] =  15.0
        x_guess[4] = -25.0
        x_guess[5] =  30.0

        tol = 5e-11
        max_it = 20

        self.coupling_utility.Initialize()

        nl_it = 1
        convergence = False

        self.coupling_utility.InitializeSolutionStep()

        while (nl_it <= max_it):

            residual = self.ComputeMVQNResidual(x_guess)
            res_norm = 0.0
            for i in range(0,len(residual)):
                res_norm += residual[i]**2
            res_norm = math.sqrt(res_norm)

            print("Iteration: ", nl_it," residual norm: ", res_norm)

            if res_norm > tol:
                self.coupling_utility.InitializeNonLinearIteration()
                self.coupling_utility.UpdateSolution(residual, x_guess)
                self.coupling_utility.FinalizeNonLinearIteration()
                nl_it += 1
            else:
                self.coupling_utility.FinalizeSolutionStep()
                convergence = True
                break

        x_final = x_guess

        # Check the obtained solution
        expected_RHS = KratosMultiphysics.Vector(6)
        expected_RHS[0] = -3.49458887
        expected_RHS[1] = 11.98381557
        expected_RHS[2] = -9.42225393
        expected_RHS[3] = 3.79688473
        expected_RHS[4] = 0.72338935
        expected_RHS[5] = 4.36818068

        obtained_RHS = self.ComputeMVQNObtainededRHS(x_final)

        if convergence == True:
            for i in range(0,len(expected_RHS)):
                self.assertAlmostEqual(expected_RHS[i],obtained_RHS[i])

    # MVQN randomized SVD accelerator test
    def testMVQNRandomizedSVD(self):

        convergence_accelerator_settings = KratosMultiphysics.Parameters("""{
            "solver_type" : "MVQN_randomized_SVD",
            "jacobian_modes" : 3,
            "min_rand_svd_extra_modes" : 2,
            "automatic_jacobian_modes" : false
        }""")

        # Construct the accelerator strategy
        self.coupling_utility = convergence_accelerator_factory.CreateConvergenceAccelerator(convergence_accelerator_settings)

        x_guess = KratosMultiphysics.Vector(10)
        for i in range(10):
            x_guess[i] = 0.0

        self.coupling_utility.Initialize()

        tol = 5e-11
        max_it = 20
        total_steps = 3
        for step in range(total_steps):
            nl_it = 1
            convergence = False
            self.coupling_utility.InitializeSolutionStep()
            while (nl_it <= max_it):
                residual = self.ComputeMVQNRandomizedSVDResidual(x_guess, step+1)
                res_norm = self.norm(residual)

                if res_norm > tol:
                    self.coupling_utility.InitializeNonLinearIteration()
                    self.coupling_utility.UpdateSolution(residual, x_guess)
                    self.coupling_utility.FinalizeNonLinearIteration()
                    nl_it += 1
                else:
                    self.coupling_utility.FinalizeSolutionStep()
                    convergence = True
                    break

        # Check the obtained solution
        expected_solution = [-0.274582199767, 0.345871967131, -3.755454999722, 2.709011159847, -4.425545499972, 1.364499472545, 8.220809505302, -1.110404752651, 9.165482205319, -9.330964410638]
        self.assertTrue(convergence)
        self.assertVectorAlmostEqual(expected_solution, x_guess)

    # Recursive MVQN accelerator test
    def test_mvqn_recusive_accelerator(self):

        mvqn_recursive_settings = KratosMultiphysics.Parameters("""{
                                                                    "solver_type" : "MVQN_recursive",
                                                                    "w_0"         : 0.825,
                                                                    "buffer_size" : 7
                                                                   }""")

        settings_list = [mvqn_recursive_settings]

        step = 1
        tol = 1e-14
        max_it = 20
        recursive_steps = 10

        x_guess = KratosMultiphysics.Vector(6)
        x_guess[0] = -20.0 * step/recursive_steps
        x_guess[1] =  20.0 * step/recursive_steps
        x_guess[2] =  -5.0 * step/recursive_steps
        x_guess[3] =  15.0 * step/recursive_steps
        x_guess[4] = -25.0 * step/recursive_steps
        x_guess[5] =  30.0 * step/recursive_steps

        while (step <= recursive_steps):

            # Iterate throug the settings list
            for i in range(0,len(settings_list)):
                print("")
                print("Testing recursive accelerator: ",settings_list[i]["solver_type"].GetString()," Step: ",step)

                # Construct the accelerator strategy
                self.coupling_utility = convergence_accelerator_factory.CreateConvergenceAccelerator(settings_list[i])
                self.coupling_utility.Initialize()

                nl_it = 1
                res_norm = 1.0
                convergence = False

                self.coupling_utility.InitializeSolutionStep()

                while (nl_it <= max_it):

                    residual = self.ComputeMVQNResidual(x_guess, step/recursive_steps)
                    res_norm = math.sqrt(residual[0]**2 + residual[1]**2)

                    print("Iteration: ", nl_it," residual norm: ", res_norm)

                    if res_norm > tol:
                        self.coupling_utility.InitializeNonLinearIteration()
                        self.coupling_utility.UpdateSolution(residual, x_guess)
                        self.coupling_utility.FinalizeNonLinearIteration()
                        nl_it += 1
                    else:
                        self.coupling_utility.FinalizeSolutionStep()
                        convergence = True
                        break

                x_final = x_guess

                # Check the obtained solution
                expected_RHS = KratosMultiphysics.Vector(6)
                expected_RHS[0] = -3.49458887 * step/recursive_steps
                expected_RHS[1] = 11.98381557 * step/recursive_steps
                expected_RHS[2] = -9.42225393 * step/recursive_steps
                expected_RHS[3] = 3.79688473 * step/recursive_steps
                expected_RHS[4] = 0.72338935 * step/recursive_steps
                expected_RHS[5] = 4.36818068 * step/recursive_steps

                obtained_RHS = self.ComputeMVQNObtainededRHS(x_final)

                if convergence == True:
                    for i in range(0,len(expected_RHS)):
                        self.assertAlmostEqual(expected_RHS[i],obtained_RHS[i])

                step += 1

    def test_accelerator_with_jacobian(self):

        settings = KratosMultiphysics.Parameters("""{
                                                     "solver_type" : "MVQN_recursive",
                                                     "w_0"         : 0.825
                                                    }""")
        print("")
        print("Testing Newton-Raphson with accelerator: ",settings["solver_type"].GetString())

        # Construct the coupling partitioned strategy
        self.coupling_utility = convergence_accelerator_factory.CreateConvergenceAccelerator(settings)

        x_guess = KratosMultiphysics.Vector(2)
        x_guess[0] = 0.5
        x_guess[1] = 0.5

        dx = KratosMultiphysics.Vector(2)
        dx[0] = 0.0
        dx[1] = 0.0

        tol = 1e-14
        max_it = 20
        acceleration = True

        if acceleration == True:
            self.coupling_utility.Initialize()

        nl_it = 1
        res_norm = 1.0
        convergence = False

        if acceleration == True:
            self.coupling_utility.InitializeSolutionStep()

        while (nl_it <= max_it):

            residual = self.ComputeResidual(x_guess)
            res_norm = math.sqrt(residual[0]**2 + residual[1]**2)

            print("Iteration: ", nl_it," residual norm: ", res_norm)

            K = self.Jacobian(x_guess)
            Kinv = self.InvertMatrix(K)
            dx = self.prod(Kinv,residual)

            if res_norm > tol:
                x_guess = x_guess + dx

                if acceleration == True:
                    residual = self.ComputeResidual(x_guess)

                    self.coupling_utility.InitializeNonLinearIteration()
                    self.coupling_utility.UpdateSolution(residual, x_guess)
                    self.coupling_utility.FinalizeNonLinearIteration()

                    residual = self.ComputeResidual(x_guess)
                    res_norm = math.sqrt( residual[0]**2 + residual[1]**2)
                    print("Iteration: ", nl_it," residual norm after acceleration: ", res_norm)

                    if res_norm < tol:
                        if acceleration == True:
                            self.coupling_utility.FinalizeSolutionStep()
                        convergence = True
                        break

                nl_it += 1
            else:
                if acceleration == True:
                    self.coupling_utility.FinalizeSolutionStep()
                convergence = True
                break

        x_final = x_guess

        # Check the obtained solution
        expected_x = KratosMultiphysics.Vector(2)
        expected_x[0] = math.sqrt(2)/2.0
        expected_x[1] = math.sqrt(2)/2.0

        if convergence == True:
            for i in range(0,len(expected_x)):
                self.assertAlmostEqual(expected_x[i],x_final[i])


    def test_accelerator_dx(self):

        settings = KratosMultiphysics.Parameters("""{
                                                     "solver_type" : "MVQN_recursive",
                                                     "w_0"         : 0.825
                                                    }""")
        print("")
        print("Testing Newton-Raphson with accelerator: ",settings["solver_type"].GetString(), " version with residual=DX")

        # Construct the coupling partitioned strategy
        self.coupling_utility = convergence_accelerator_factory.CreateConvergenceAccelerator(settings)

        x_guess = KratosMultiphysics.Vector(2)
        x_guess[0] = 0.5
        x_guess[1] = 0.5

        dx = KratosMultiphysics.Vector(2)
        dx[0] = 0.0
        dx[1] = 0.0

        tol = 1e-14
        max_it = 20
        acceleration = True

        if acceleration == True:
            self.coupling_utility.Initialize()

        nl_it = 1
        res_norm = 1.0
        convergence = False

        if acceleration == True:
            self.coupling_utility.InitializeSolutionStep()

        while (nl_it <= max_it):

            residual = self.ComputeResidual(x_guess)
            res_norm = math.sqrt(residual[0]**2 + residual[1]**2)

            print("Iteration: ", nl_it," residual norm: ", res_norm)

            K = self.Jacobian(x_guess)
            Kinv = self.InvertMatrix(K)
            dx = self.prod(Kinv,residual)

            if res_norm > tol:
                #x_guess = x_guess + dx

                if acceleration == True:
                    residual = dx #self.ComputeResidual(x_guess)

                    self.coupling_utility.InitializeNonLinearIteration()
                    self.coupling_utility.UpdateSolution(residual, x_guess)
                    self.coupling_utility.FinalizeNonLinearIteration()

                    residual = self.ComputeResidual(x_guess)
                    res_norm = math.sqrt( residual[0]**2 + residual[1]**2)
                    print("Iteration: ", nl_it," residual norm after acceleration: ", res_norm)

                    if res_norm < tol:
                        if acceleration == True:
                            self.coupling_utility.FinalizeSolutionStep()
                        convergence = True
                        break

                nl_it += 1
            else:
                if acceleration == True:
                    self.coupling_utility.FinalizeSolutionStep()
                convergence = True
                break

        x_final = x_guess

        # Check the obtained solution
        expected_x = KratosMultiphysics.Vector(2)
        expected_x[0] = math.sqrt(2)/2.0
        expected_x[1] = math.sqrt(2)/2.0

        if convergence == True:
            for i in range(0,len(expected_x)):
                self.assertAlmostEqual(expected_x[i],x_final[i])

if __name__ == '__main__':
    KratosUnittest.main()
