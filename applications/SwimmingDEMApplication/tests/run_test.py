# import Kratos
import sys

import KratosMultiphysics
import KratosMultiphysics.ExternalSolversApplication
import KratosMultiphysics.FluidDynamicsApplication
import KratosMultiphysics.DEMApplication
import KratosMultiphysics.SwimmingDEMApplication

# Import Kratos "wrapper" for unittests
import KratosMultiphysics.KratosUnittest as KratosUnittest

# Import the tests o test_classes to create the suits:

# SMALL TESTS
import SmallTests
# NIGTHLY TESTS
import NightTests

if __name__ == '__main__':
    #case = eval(sys.argv[1])
    case = test_SwimmingDEMApplication

    KratosUnittest.TextTestRunner(verbosity=2).run(
        KratosUnittest.TestLoader().loadTestsFromTestCase(case)
    )
