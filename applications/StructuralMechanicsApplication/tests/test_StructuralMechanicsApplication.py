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
# from SmallTests import SimpleMeshMovingTest             as TSimpleMeshMovingTest
# # Dynamic basic tests
# from SmallTests import DynamicBossakTests               as TDynamicBossakTests
# from SmallTests import DynamicNewmarkTests              as TDynamicNewmarkTests
# # Patch test Small Displacements
# from SmallTests import SDTwoDShearQuaPatchTest          as TSDTwoDShearQuaPatchTest
# from SmallTests import SDTwoDShearTriPatchTest          as TSDTwoDShearTriPatchTest
# from SmallTests import SDTwoDTensionQuaPatchTest        as TSDTwoDTensionQuaPatchTest
# from SmallTests import SDTwoDTensionTriPatchTest        as TSDTwoDTensionTriPatchTest
# from SmallTests import SDThreeDShearHexaPatchTest       as TSDThreeDShearHexaPatchTest
# from SmallTests import SDThreeDShearTetraPatchTest      as TSDThreeDShearTetraPatchTest
# from SmallTests import SDThreeDTensionHexaPatchTest     as TSDThreeDTensionHexaPatchTest
# from SmallTests import SDThreeDTensionTetraPatchTest    as TSDThreeDTensionTetraPatchTest
# # Patch test Total Lagrangian
# from SmallTests import TLTwoDShearQuaPatchTest          as TTLTwoDShearQuaPatchTest
# from SmallTests import TLTwoDShearTriPatchTest          as TTLTwoDShearTriPatchTest
# from SmallTests import TLTwoDTensionQuaPatchTest        as TTLTwoDTensionQuaPatchTest
# from SmallTests import TLTwoDTensionTriPatchTest        as TTLTwoDTensionTriPatchTest
# from SmallTests import TLThreeDShearHexaPatchTest       as TTLThreeDShearHexaPatchTest
# from SmallTests import TLThreeDShearTetraPatchTest      as TTLThreeDShearTetraPatchTest
# from SmallTests import TLThreeDTensionHexaPatchTest     as TTLThreeDTensionHexaPatchTest
# from SmallTests import TLThreeDTensionTetraPatchTest    as TTLThreeDTensionTetraPatchTest
# # Patch test Updated Lagrangian
# from SmallTests import ULTwoDShearQuaPatchTest          as TULTwoDShearQuaPatchTest
# from SmallTests import ULTwoDShearTriPatchTest          as TULTwoDShearTriPatchTest
# from SmallTests import ULTwoDTensionQuaPatchTest        as TULTwoDTensionQuaPatchTest
# from SmallTests import ULTwoDTensionTriPatchTest        as TULTwoDTensionTriPatchTest
# from SmallTests import ULThreeDShearHexaPatchTest       as TULThreeDShearHexaPatchTest
# from SmallTests import ULThreeDShearTetraPatchTest      as TULThreeDShearTetraPatchTest
# from SmallTests import ULThreeDTensionHexaPatchTest     as TULThreeDTensionHexaPatchTest
# from SmallTests import ULThreeDTensionTetraPatchTest    as TULThreeDTensionTetraPatchTest
# # SPRISM tests
# from SmallTests import SprismMembranePatchTests         as TSprismMembranePatchTests
# from SmallTests import SprismBendingPatchTests          as TSprismBendingPatchTests
# # Shell tests
# from SmallTests import ShellQ4ThickBendingRollUpTests   as TShellQ4ThickBendingRollUpTests
# from SmallTests import ShellQ4ThickDrillingRollUpTests  as TShellQ4ThickDrillingRollUpTests
from SmallTests import ShellT3ThinBendingRollUpTests    as TShellT3ThinBendingRollUpTests
# from SmallTests import ShellT3ThinDrillingRollUpTests   as TShellT3ThinDrillingRollUpTests
# # Eigenvalues tests
# from SmallTests import EigenQ4Thick2x2PlateTests        as TEigenQ4Thick2x2PlateTests
# from SmallTests import EigenTL3D8NCubeTests             as TEigenTL3D8NCubeTests
# # Membrane tests
# from SmallTests import Fofi4PointTentnoCableTests       as TFofi4PointTentnoCableTests
from SmallTests import MembraneQ4PointLoadTests         as TMembraneQ4PointLoadTests
# # Nodal damping test
# from test_nodal_damping import NodalDampingTests        as TNodalDampingTests
# 2NELEMENT tests
from SmallTests import Simple3D2NTrussTest as T3D2NTrussTest
from SmallTests import Simple3D2NTrussLinearTest as T3D2NTrussLinearTest
from SmallTests import Simple3D2NTrussDynamicTest as T3D2NTrussDynamicTest
from SmallTests import Simple3D2NBeamCrTest as T3D2NBeamCrTest
from SmallTests import Simple3D2NBeamCrLinearTest as T3D2NBeamCrLinearTest
from SmallTests import Simple3D2NBeamCrDynamicTest as T3D2NBeamCrDynamicTest
# Nodal damping test
from test_nodal_damping import NodalDampingTests        as TNodalDampingTests
from test_spring_damper_element import SpringDamperElementTests as TSpringDamperElementTests

# ## NIGTHLY TESTS
# from NightlyTests import ShellT3IsotropicScordelisTests as TShellT3IsotropicScordelisTests

# ## VALIDATION TESTS
# from ValidationTests import SprismPanTests as TSprismPanTests

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
    # smallSuite.addTest(TSimpleMeshMovingTest('test_execution'))
    # smallSuite.addTest(TDynamicBossakTests('test_execution'))
    # smallSuite.addTest(TDynamicNewmarkTests('test_execution'))
    # smallSuite.addTest(TSprismMembranePatchTests('test_execution'))
    # smallSuite.addTest(TSprismBendingPatchTests('test_execution'))
    # smallSuite.addTest(TShellQ4ThickBendingRollUpTests('test_execution'))
    # smallSuite.addTest(TShellQ4ThickDrillingRollUpTests('test_execution'))
    smallSuite.addTest(TShellT3ThinBendingRollUpTests('test_execution'))
    # smallSuite.addTest(TShellT3ThinDrillingRollUpTests('test_execution'))
    # smallSuite.addTest(TEigenQ4Thick2x2PlateTests('test_execution'))
    # smallSuite.addTest(TEigenTL3D8NCubeTests('test_execution'))
    
    smallSuite.addTest(T3D2NTrussDynamicTest('test_execution')) 
    smallSuite.addTest(T3D2NTrussLinearTest('test_execution'))  
    smallSuite.addTest(T3D2NTrussTest('test_execution'))  
    smallSuite.addTest(T3D2NBeamCrTest('test_execution'))    
    smallSuite.addTest(T3D2NBeamCrLinearTest('test_execution'))  
    smallSuite.addTest(T3D2NBeamCrDynamicTest('test_execution'))  

    # smallSuite.addTest(MPCSmallDisplacementElementTests('test_execution'))

    # # Membrane tests
    # smallSuite.addTest(TFofi4PointTentnoCableTests('test_execution'))
    smallSuite.addTest(TMembraneQ4PointLoadTests('test_execution'))
    # # Nodal damping test
    # smallSuite.addTest(TNodalDampingTests('test_execution'))
    # Nodal damping test
    smallSuite.addTest(TNodalDampingTests('test_execution'))
    smallSuite.addTest(TSpringDamperElementTests('test_execution'))

    # Create a test suit with the selected tests plus all small tests
    nightSuite = suites['nightly']
    nightSuite.addTests(smallSuite)
    #nightSuite.addTest(TShellT3IsotropicScordelisTests('test_execution'))
    
    # For very long tests that should not be in nighly and you can use to validate 
    validationSuite = suites['validation']
    #validationSuite.addTest(TSprismPanTests('test_execution'))

    # Create a test suit that contains all the tests:
    allSuite = suites['all']
    allSuite.addTests(
        KratosUnittest.TestLoader().loadTestsFromTestCases([
            #TFofi4PointTentnoCableTests,
            #TMembraneQ4PointLoadTests,
            #TSimpleMeshMovingTest,
            #TDynamicBossakTests,
            #TDynamicNewmarkTests,
            #TSprismMembranePatchTests,
            #TSprismBendingPatchTests,
            #TShellQ4ThickBendingRollUpTests,
            #TShellQ4ThickDrillingRollUpTests,
            #TShellT3ThinBendingRollUpTests,
            #TShellT3ThinDrillingRollUpTests,
            #TShellT3IsotropicScordelisTests,
            #TMembraneQ4PointLoadTests,
            T3D2NTrussTest,
            T3D2NTrussLinearTest,
            T3D2NTrussDynamicTest,
            T3D2NBeamCrTest,
            T3D2NBeamCrLinearTest,
            T3D2NBeamCrDynamicTest
            #TMpcSmallDispElemTests
            #TIsotropicDamageSimoJuPSTest,
            #TNodalDampingTests
            #TShellT3ThinDrillingRollUpTests,
            #TShellT3IsotropicScordelisTests,
            #TIsotropicDamageSimoJuPSTest,
            #TNodalDampingTests,
            #TSpringDamperElementTests
            ######TSprismPanTests
        ])
    )
    
    if( hasattr(KratosMultiphysics.ExternalSolversApplication,  "FEASTSolver") ):
        allSuite.addTests(
            KratosUnittest.TestLoader().loadTestsFromTestCases([
               # TEigenQ4Thick2x2PlateTests,
               # TEigenTL3D8NCubeTests
            ])
        )
    else:
        print("FEASTSolver solver is not included in the compilation of the External Solvers Application")

    return suites

if __name__ == '__main__':
    KratosUnittest.runTests(AssambleTestSuites())
