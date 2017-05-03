# import Kratos
import KratosMultiphysics 
import KratosMultiphysics.ExternalSolversApplication 
import KratosMultiphysics.SolidMechanicsApplication 

# Import Kratos "wrapper" for unittests
import KratosMultiphysics.KratosUnittest as KratosUnittest

# Import the tests o test_classes to create the suits
## SMALL TESTS
# Patch test Small Displacements
from SmallTests import SDTwoDShearQuaPatchTest       as TSDTwoDShearQuaPatchTest
from SmallTests import SDTwoDShearTriPatchTest       as TSDTwoDShearTriPatchTest
from SmallTests import SDTwoDTensionQuaPatchTest     as TSDTwoDTensionQuaPatchTest
from SmallTests import SDTwoDTensionTriPatchTest     as TSDTwoDTensionTriPatchTest
from SmallTests import SDThreeDShearHexaPatchTest    as TSDThreeDShearHexaPatchTest
from SmallTests import SDThreeDShearTetraPatchTest   as TSDThreeDShearTetraPatchTest
from SmallTests import SDThreeDTensionHexaPatchTest  as TSDThreeDTensionHexaPatchTest
from SmallTests import SDThreeDTensionTetraPatchTest as TSDThreeDTensionTetraPatchTest
# Patch test Total Lagrangian
from SmallTests import TLTwoDShearQuaPatchTest       as TTLTwoDShearQuaPatchTest
from SmallTests import TLTwoDShearTriPatchTest       as TTLTwoDShearTriPatchTest
from SmallTests import TLTwoDTensionQuaPatchTest     as TTLTwoDTensionQuaPatchTest
from SmallTests import TLTwoDTensionTriPatchTest     as TTLTwoDTensionTriPatchTest
from SmallTests import TLThreeDShearHexaPatchTest    as TTLThreeDShearHexaPatchTest
from SmallTests import TLThreeDShearTetraPatchTest   as TTLThreeDShearTetraPatchTest
from SmallTests import TLThreeDTensionHexaPatchTest  as TTLThreeDTensionHexaPatchTest
from SmallTests import TLThreeDTensionTetraPatchTest as TTLThreeDTensionTetraPatchTest
# Patch test Updated Lagrangian
from SmallTests import ULTwoDShearQuaPatchTest       as TULTwoDShearQuaPatchTest
from SmallTests import ULTwoDShearTriPatchTest       as TULTwoDShearTriPatchTest
from SmallTests import ULTwoDTensionQuaPatchTest     as TULTwoDTensionQuaPatchTest
from SmallTests import ULTwoDTensionTriPatchTest     as TULTwoDTensionTriPatchTest
from SmallTests import ULThreeDShearHexaPatchTest    as TULThreeDShearHexaPatchTest
from SmallTests import ULThreeDShearTetraPatchTest   as TULThreeDShearTetraPatchTest
from SmallTests import ULThreeDTensionHexaPatchTest  as TULThreeDTensionHexaPatchTest
from SmallTests import ULThreeDTensionTetraPatchTest as TULThreeDTensionTetraPatchTest

## NIGTHLY TESTS

## VALIDATION TESTS

def AssambleTestSuites():
    ''' Populates the test suites to run.

    Populates the test suites to run. At least, it should pupulate the suites:
    "small", "nighlty" and "all"

    Return
    ------

    suites: A dictionary of suites
        The set of suites with its test_cases added.
    '''
    suites = KratosUnittest.KratosSuites

    # Create a test suit with the selected tests (Small tests):
    smallSuite = suites['small']
    smallSuite.addTest(TSDTwoDShearQuaPatchTest('test_execution'))
    smallSuite.addTest(TSDTwoDShearTriPatchTest('test_execution'))
    smallSuite.addTest(TSDTwoDTensionQuaPatchTest('test_execution'))
    smallSuite.addTest(TSDTwoDTensionTriPatchTest('test_execution'))
    smallSuite.addTest(TSDThreeDShearHexaPatchTest('test_execution'))
    smallSuite.addTest(TSDThreeDShearTetraPatchTest('test_execution'))
    smallSuite.addTest(TSDThreeDTensionHexaPatchTest('test_execution'))
    smallSuite.addTest(TSDThreeDTensionTetraPatchTest('test_execution'))
    smallSuite.addTest(TTLTwoDShearQuaPatchTest('test_execution'))
    smallSuite.addTest(TTLTwoDShearTriPatchTest('test_execution'))
    smallSuite.addTest(TTLTwoDTensionQuaPatchTest('test_execution'))
    smallSuite.addTest(TTLTwoDTensionTriPatchTest('test_execution'))
    smallSuite.addTest(TTLThreeDShearHexaPatchTest('test_execution'))
    smallSuite.addTest(TTLThreeDShearTetraPatchTest('test_execution'))
    smallSuite.addTest(TTLThreeDTensionHexaPatchTest('test_execution'))
    smallSuite.addTest(TTLThreeDTensionTetraPatchTest('test_execution'))
    smallSuite.addTest(TULTwoDShearQuaPatchTest('test_execution'))
    smallSuite.addTest(TULTwoDShearTriPatchTest('test_execution'))
    smallSuite.addTest(TULTwoDTensionQuaPatchTest('test_execution'))
    smallSuite.addTest(TULTwoDTensionTriPatchTest('test_execution'))
    smallSuite.addTest(TULThreeDShearHexaPatchTest('test_execution'))
    smallSuite.addTest(TULThreeDShearTetraPatchTest('test_execution'))
    smallSuite.addTest(TULThreeDTensionHexaPatchTest('test_execution'))
    smallSuite.addTest(TULThreeDTensionTetraPatchTest('test_execution'))
    
    # Create a test suit with the selected tests plus all small tests
    nightSuite = suites['nightly']
    nightSuite.addTests(smallSuite)
    
    # For very long tests that should not be in nighly and you can use to validate 
    validationSuite = suites['validation']

    # Create a test suit that contains all the tests:
    allSuite = suites['all']
    allSuite.addTests(
        KratosUnittest.TestLoader().loadTestsFromTestCases([
            TSDTwoDShearQuaPatchTest,
            TSDTwoDShearTriPatchTest,
            TSDTwoDTensionQuaPatchTest,
            TSDTwoDTensionTriPatchTest,
            TSDThreeDShearHexaPatchTest,
            TSDThreeDShearTetraPatchTest,
            TSDThreeDTensionHexaPatchTest,
            TSDThreeDTensionTetraPatchTest,
            TTLTwoDShearQuaPatchTest,
            TTLTwoDShearTriPatchTest,
            TTLTwoDTensionQuaPatchTest,
            TTLTwoDTensionTriPatchTest,
            TTLThreeDShearHexaPatchTest,
            TTLThreeDShearTetraPatchTest,
            TTLThreeDTensionHexaPatchTest,
            TTLThreeDTensionTetraPatchTest,
            TULTwoDShearQuaPatchTest,
            TULTwoDShearTriPatchTest,
            TULTwoDTensionQuaPatchTest,
            TULTwoDTensionTriPatchTest,
            TULThreeDShearHexaPatchTest,
            TULThreeDShearTetraPatchTest,
            TULThreeDTensionHexaPatchTest,
            TULThreeDTensionTetraPatchTest,
        ])
    )

    return suites

if __name__ == '__main__':
    KratosUnittest.runTests(AssambleTestSuites())
