# import Kratos
import sys

import KratosMultiphysics
import KratosMultiphysics.ConstitutiveModelsApplication

# Import Kratos "wrapper" for unittests
import KratosMultiphysics.KratosUnittest as KratosUnittest

# Import the tests o test_classes to create the suits:

# VALIDATION TESTS
import ValidationTests


if __name__ == '__main__':
    case = eval(sys.argv[1])

    KratosUnittest.TextTestRunner(verbosity=2).run(
        KratosUnittest.TestLoader().loadTestsFromTestCase(case)
    )
