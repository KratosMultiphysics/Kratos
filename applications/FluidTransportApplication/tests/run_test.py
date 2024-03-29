# import Kratos
import sys

import KratosMultiphysics
import KratosMultiphysics.LinearSolversApplication
import KratosMultiphysics.ConvectionDiffusionApplication
import KratosMultiphysics.FluidDynamicsApplication
import KratosMultiphysics.FluidTransportApplication

# Import Kratos "wrapper" for unittests
import KratosMultiphysics.KratosUnittest as KratosUnittest

# Import the tests o test_classes to create the suits:

# SMALL TESTS
import SmallTests
# NIGTHLY TESTS
import NightTests

if __name__ == '__main__':
    case = eval(sys.argv[1])

    KratosUnittest.TextTestRunner(verbosity=2).run(
        KratosUnittest.TestLoader().loadTestsFromTestCase(case)
    )
