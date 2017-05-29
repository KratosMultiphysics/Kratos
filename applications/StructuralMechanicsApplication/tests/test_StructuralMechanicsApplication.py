# import Kratos
import KratosMultiphysics 
import KratosMultiphysics.ExternalSolversApplication 
import KratosMultiphysics.SolidMechanicsApplication 
import KratosMultiphysics.StructuralMechanicsApplication 

# Import Kratos "wrapper" for unittests
import KratosMultiphysics.KratosUnittest as KratosUnittest

# Import the tests o test_classes to create the suits
## SMALL TESTS

# Basic moving mesh test
from SmallTests import SimpleMeshMovingTest             as TSimpleMeshMovingTest
# Dynamic basic tests
from SmallTests import DynamicBossakTests               as TDynamicBossakTests
from SmallTests import DynamicNewmarkTests              as TDynamicNewmarkTests
# Patch test Small Displacements
from SmallTests import SDTwoDShearQuaPatchTest          as TSDTwoDShearQuaPatchTest
from SmallTests import SDTwoDShearTriPatchTest          as TSDTwoDShearTriPatchTest
from SmallTests import SDTwoDTensionQuaPatchTest        as TSDTwoDTensionQuaPatchTest
from SmallTests import SDTwoDTensionTriPatchTest        as TSDTwoDTensionTriPatchTest
from SmallTests import SDThreeDShearHexaPatchTest       as TSDThreeDShearHexaPatchTest
from SmallTests import SDThreeDShearTetraPatchTest      as TSDThreeDShearTetraPatchTest
from SmallTests import SDThreeDTensionHexaPatchTest     as TSDThreeDTensionHexaPatchTest
from SmallTests import SDThreeDTensionTetraPatchTest    as TSDThreeDTensionTetraPatchTest
# Patch test Total Lagrangian
from SmallTests import TLTwoDShearQuaPatchTest          as TTLTwoDShearQuaPatchTest
from SmallTests import TLTwoDShearTriPatchTest          as TTLTwoDShearTriPatchTest
from SmallTests import TLTwoDTensionQuaPatchTest        as TTLTwoDTensionQuaPatchTest
from SmallTests import TLTwoDTensionTriPatchTest        as TTLTwoDTensionTriPatchTest
from SmallTests import TLThreeDShearHexaPatchTest       as TTLThreeDShearHexaPatchTest
from SmallTests import TLThreeDShearTetraPatchTest      as TTLThreeDShearTetraPatchTest
from SmallTests import TLThreeDTensionHexaPatchTest     as TTLThreeDTensionHexaPatchTest
from SmallTests import TLThreeDTensionTetraPatchTest    as TTLThreeDTensionTetraPatchTest
# Patch test Updated Lagrangian
from SmallTests import ULTwoDShearQuaPatchTest          as TULTwoDShearQuaPatchTest
from SmallTests import ULTwoDShearTriPatchTest          as TULTwoDShearTriPatchTest
from SmallTests import ULTwoDTensionQuaPatchTest        as TULTwoDTensionQuaPatchTest
from SmallTests import ULTwoDTensionTriPatchTest        as TULTwoDTensionTriPatchTest
from SmallTests import ULThreeDShearHexaPatchTest       as TULThreeDShearHexaPatchTest
from SmallTests import ULThreeDShearTetraPatchTest      as TULThreeDShearTetraPatchTest
from SmallTests import ULThreeDTensionHexaPatchTest     as TULThreeDTensionHexaPatchTest
from SmallTests import ULThreeDTensionTetraPatchTest    as TULThreeDTensionTetraPatchTest
# SPRISM tests
from SmallTests import SprismMembranePatchTests         as TSprismMembranePatchTests
from SmallTests import SprismBendingPatchTests          as TSprismBendingPatchTests
# Shell tests
from SmallTests import ShellQ4ThickBendingRollUpTests   as TShellQ4ThickBendingRollUpTests
from SmallTests import ShellQ4ThickDrillingRollUpTests  as TShellQ4ThickDrillingRollUpTests
from SmallTests import ShellT3ThinBendingRollUpTests    as TShellT3ThinBendingRollUpTests
from SmallTests import ShellT3ThinDrillingRollUpTests   as TShellT3ThinDrillingRollUpTests
# Eigenvalues tests
from SmallTests import EigenQ4Thick2x2PlateTests        as TEigenQ4Thick2x2PlateTests
from SmallTests import EigenTL3D8NCubeTests             as TEigenTL3D8NCubeTests
# Membrane tests
from SmallTests import Fofi4PointTentnoCableTests       as TFofi4PointTentnoCableTests
from SmallTests import MembraneQ4PointLoadTests         as TMembraneQ4PointLoadTests
# Nodal damping test
from test_nodal_damping import NodalDampingTests        as TNodalDampingTests

## NIGTHLY TESTS
# Shell test
from NightlyTests import ShellT3IsotropicScordelisTests as TShellT3IsotropicScordelisTests
# CL tests
from NightlyTests import IsotropicDamageSimoJuPSTest    as TIsotropicDamageSimoJuPSTest

## VALIDATION TESTS
# SPRISM tests
from ValidationTests import SprismPanTests              as TSprismPanTests
# Eigenvalues tests
from ValidationTests import Eigen3D3NThinCircleTests    as TEigen3D3NThinCircleTests

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
    # Basic moving mesh test
    smallSuite.addTest(TSimpleMeshMovingTest('test_execution'))
    # Dynamic basic tests
    smallSuite.addTest(TDynamicBossakTests('test_execution'))
    smallSuite.addTest(TDynamicNewmarkTests('test_execution'))
    # Patch test Small Displacements
    smallSuite.addTest(TSDTwoDShearQuaPatchTest('test_execution'))
    smallSuite.addTest(TSDTwoDShearTriPatchTest('test_execution'))
    smallSuite.addTest(TSDTwoDTensionQuaPatchTest('test_execution'))
    smallSuite.addTest(TSDTwoDTensionTriPatchTest('test_execution'))
    smallSuite.addTest(TSDThreeDShearHexaPatchTest('test_execution'))
    smallSuite.addTest(TSDThreeDShearTetraPatchTest('test_execution'))
    smallSuite.addTest(TSDThreeDTensionHexaPatchTest('test_execution'))
    smallSuite.addTest(TSDThreeDTensionTetraPatchTest('test_execution'))
    # Patch test Total Lagrangian
    smallSuite.addTest(TTLTwoDShearQuaPatchTest('test_execution'))
    smallSuite.addTest(TTLTwoDShearTriPatchTest('test_execution'))
    smallSuite.addTest(TTLTwoDTensionQuaPatchTest('test_execution'))
    smallSuite.addTest(TTLTwoDTensionTriPatchTest('test_execution'))
    smallSuite.addTest(TTLThreeDShearHexaPatchTest('test_execution'))
    smallSuite.addTest(TTLThreeDShearTetraPatchTest('test_execution'))
    smallSuite.addTest(TTLThreeDTensionHexaPatchTest('test_execution'))
    smallSuite.addTest(TTLThreeDTensionTetraPatchTest('test_execution'))
    # Patch test Updated Lagrangian
    smallSuite.addTest(TULTwoDShearQuaPatchTest('test_execution'))
    smallSuite.addTest(TULTwoDShearTriPatchTest('test_execution'))
    smallSuite.addTest(TULTwoDTensionQuaPatchTest('test_execution'))
    smallSuite.addTest(TULTwoDTensionTriPatchTest('test_execution'))
    smallSuite.addTest(TULThreeDShearHexaPatchTest('test_execution'))
    smallSuite.addTest(TULThreeDShearTetraPatchTest('test_execution'))
    smallSuite.addTest(TULThreeDTensionHexaPatchTest('test_execution'))
    smallSuite.addTest(TULThreeDTensionTetraPatchTest('test_execution'))
    # SPRISM tests
    smallSuite.addTest(TSprismMembranePatchTests('test_execution'))
    smallSuite.addTest(TSprismBendingPatchTests('test_execution'))
    # Shell tests
    smallSuite.addTest(TShellQ4ThickBendingRollUpTests('test_execution'))
    smallSuite.addTest(TShellQ4ThickDrillingRollUpTests('test_execution'))
    smallSuite.addTest(TShellT3ThinBendingRollUpTests('test_execution'))
    smallSuite.addTest(TShellT3ThinDrillingRollUpTests('test_execution'))
    # Eigenvalues tests
    smallSuite.addTest(TEigenQ4Thick2x2PlateTests('test_execution'))
    smallSuite.addTest(TEigenTL3D8NCubeTests('test_execution'))
    # Membrane tests
    smallSuite.addTest(TFofi4PointTentnoCableTests('test_execution'))
    smallSuite.addTest(TMembraneQ4PointLoadTests('test_execution'))
    # Nodal damping test
    smallSuite.addTest(TNodalDampingTests('test_execution'))

    # Create a test suit with the selected tests plus all small tests
    nightSuite = suites['nightly']
    nightSuite.addTests(smallSuite)
    # Shell tests
    nightSuite.addTest(TShellT3IsotropicScordelisTests('test_execution'))
    # CL tests
    nightSuite.addTest(TIsotropicDamageSimoJuPSTest('test_execution'))
    
    # For very long tests that should not be in nighly and you can use to validate 
    validationSuite = suites['validation']
    # SPRISM tests
    validationSuite.addTest(TSprismPanTests('test_execution'))
    # Eigenvalues tests
    validationSuite.addTest(TEigen3D3NThinCircleTests('test_execution'))

    # Create a test suit that contains all the tests:
    allSuite = suites['all']
    allSuite.addTests(
        KratosUnittest.TestLoader().loadTestsFromTestCases([
            #TFofi4PointTentnoCableTests,
            TMembraneQ4PointLoadTests,
            TSimpleMeshMovingTest,
            TDynamicBossakTests,
            TDynamicNewmarkTests,
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
            TSprismMembranePatchTests,
            TSprismBendingPatchTests,
            TShellQ4ThickBendingRollUpTests,
            TShellQ4ThickDrillingRollUpTests,
            TShellT3ThinBendingRollUpTests,
            TShellT3ThinDrillingRollUpTests,
            TShellT3IsotropicScordelisTests,
            TIsotropicDamageSimoJuPSTest,
            TNodalDampingTests
            ######TSprismPanTests
        ])
    )
    
    if( hasattr(KratosMultiphysics.ExternalSolversApplication,  "FEASTSolver") ):
        allSuite.addTests(
            KratosUnittest.TestLoader().loadTestsFromTestCases([
                TEigenQ4Thick2x2PlateTests,
                TEigenTL3D8NCubeTests,
                TEigen3D3NThinCircleTests
            ])
        )
    else:
        print("FEASTSolver solver is not included in the compilation of the External Solvers Application")

    return suites

if __name__ == '__main__':
    KratosUnittest.runTests(AssambleTestSuites())
