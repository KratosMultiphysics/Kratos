# import Kratos
import KratosMultiphysics
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication

# Import Kratos "wrapper" for unittests
import KratosMultiphysics.KratosUnittest as KratosUnittest

try:
    import KratosMultiphysics.ExternalSolversApplication as ExternalSolversApplication
    missing_external_dependencies = False
    missing_application = ''
except ImportError as e:
    missing_external_dependencies = True
    # extract name of the missing application from the error message
    import re
    missing_application = re.search(r'''.*'KratosMultiphysics\.(.*)'.*''',
                                    '{0}'.format(e)).group(1)

# Import the tests or test_classes to create the suits

##### SELF-CONTAINED TESTS #####
# CL tests
from constitutive_law_test import TestConstitutiveLaw as TTestConstitutiveLaw
# Processes test
from test_mass_calculation import TestMassCalculation as TTestMassCalculation
# Simple patch tests
from test_patch_test_small_strain import TestPatchTestSmallStrain as TTestPatchTestSmallStrain
from test_patch_test_large_strain import TestPatchTestLargeStrain as TTestPatchTestLargeStrain
from test_quadratic_elements import TestQuadraticElements as TTestQuadraticElements
from test_patch_test_shells import TestPatchTestShells as TTestPatchTestShells
from test_patch_test_truss import TestTruss3D2N as TTestTruss3D2N
from test_patch_test_cr_beam import TestCrBeam3D2N as TTestCrBeam3D2N
from test_patch_test_cr_beam import TestCrBeam2D2N as TTestCrBeam2D2N
from test_patch_test_shells_stress import TestPatchTestShellsStressRec as TTestPatchTestShellsStressRec
from test_patch_test_shells_orthotropic import TestPatchTestShellsOrthotropic as TTestPatchTestShellsOrthotropic
from test_patch_test_formfinding import TestPatchTestFormfinding as TTestPatchTestFormfinding
# Test loading conditions
from test_loading_conditions_point import TestLoadingConditionsPoint as TTestLoadingConditionsPoint
from test_loading_conditions_line import TestLoadingConditionsLine as TTestLoadingConditionsLine
from test_loading_conditions_surface import TestLoadingConditionsSurface as TTestLoadingConditionsSurface
# Multipoint constraint tests
from test_multipoint_contstraints import TestMultipointConstraints as TTestMultipointConstraints
# Nodal damping test
from test_nodal_damping import NodalDampingTests as TNodalDampingTests
# Spring damper element test
from test_spring_damper_element import SpringDamperElementTests as TSpringDamperElementTests
# Harmonic analysis tests
from test_harmonic_analysis import HarmonicAnalysisTests as THarmonicAnalysisTests


##### SMALL TESTS #####
# Dynamic basic tests (leave these in the smallSuite to have the Exection script tested)
from SmallTests import DynamicBossakTests as TDynamicBossakTests
from SmallTests import DynamicNewmarkTests as TDynamicNewmarkTests


##### NIGHTLY TESTS #####
# Basic moving mesh test
from NightlyTests import SimpleMeshMovingTest as TSimpleMeshMovingTest
# Patch test Small Displacements
from NightlyTests import SDTwoDShearQuaPatchTest as TSDTwoDShearQuaPatchTest
from NightlyTests import SDTwoDShearTriPatchTest as TSDTwoDShearTriPatchTest
from NightlyTests import SDTwoDTensionQuaPatchTest as TSDTwoDTensionQuaPatchTest
from NightlyTests import SDTwoDTensionTriPatchTest as TSDTwoDTensionTriPatchTest
from NightlyTests import SDThreeDShearHexaPatchTest as TSDThreeDShearHexaPatchTest
from NightlyTests import SDThreeDShearTetraPatchTest as TSDThreeDShearTetraPatchTest
from NightlyTests import SDThreeDTensionHexaPatchTest as TSDThreeDTensionHexaPatchTest
from NightlyTests import SDThreeDTensionTetraPatchTest as TSDThreeDTensionTetraPatchTest
# Patch test Total Lagrangian
from NightlyTests import TLTwoDShearQuaPatchTest as TTLTwoDShearQuaPatchTest
from NightlyTests import TLTwoDShearTriPatchTest as TTLTwoDShearTriPatchTest
from NightlyTests import TLTwoDTensionQuaPatchTest as TTLTwoDTensionQuaPatchTest
from NightlyTests import TLTwoDTensionTriPatchTest as TTLTwoDTensionTriPatchTest
from NightlyTests import TLThreeDShearHexaPatchTest as TTLThreeDShearHexaPatchTest
from NightlyTests import TLThreeDShearTetraPatchTest as TTLThreeDShearTetraPatchTest
from NightlyTests import TLThreeDTensionHexaPatchTest as TTLThreeDTensionHexaPatchTest
from NightlyTests import TLThreeDTensionTetraPatchTest as TTLThreeDTensionTetraPatchTest
# Patch test Updated Lagrangian
from NightlyTests import ULTwoDShearQuaPatchTest as TULTwoDShearQuaPatchTest
from NightlyTests import ULTwoDShearTriPatchTest as TULTwoDShearTriPatchTest
from NightlyTests import ULTwoDTensionQuaPatchTest as TULTwoDTensionQuaPatchTest
from NightlyTests import ULTwoDTensionTriPatchTest as TULTwoDTensionTriPatchTest
from NightlyTests import ULThreeDShearHexaPatchTest as TULThreeDShearHexaPatchTest
from NightlyTests import ULThreeDShearTetraPatchTest as TULThreeDShearTetraPatchTest
from NightlyTests import ULThreeDTensionHexaPatchTest as TULThreeDTensionHexaPatchTest
from NightlyTests import ULThreeDTensionTetraPatchTest as TULThreeDTensionTetraPatchTest
# SPRISM tests
from NightlyTests import SprismMembranePatchTests as TSprismMembranePatchTests
from NightlyTests import SprismBendingPatchTests as TSprismBendingPatchTests
# Eigenvalues tests
from NightlyTests import EigenQ4Thick2x2PlateTests as TEigenQ4Thick2x2PlateTests
from NightlyTests import EigenTL3D8NCubeTests as TEigenTL3D8NCubeTests
from NightlyTests import Eigen3D3NThinCircleTests as TEigen3D3NThinCircleTests
# Membrane tests
from NightlyTests import Fofi4PointTentnoCableTests as TFofi4PointTentnoCableTests
from NightlyTests import Fofi4PointTentCableTests as TFofi4PointTentCableTests
from NightlyTests import MembraneQ4PointLoadTests as TMembraneQ4PointLoadTests
from NightlyTests import MembraneQ4TrussPointLoadTests as TMembraneQ4TrussPointLoadTests
# 2Node Element tests
from NightlyTests import Simple3D2NTrussTest as T3D2NTrussTest
from NightlyTests import Simple3D2NTrussLinearTest as T3D2NTrussLinearTest
from NightlyTests import Simple3D2NTrussDynamicTest as T3D2NTrussDynamicTest
from NightlyTests import Simple3D2NBeamCrTest as T3D2NBeamCrTest
from NightlyTests import Simple3D2NBeamCrLinearTest as T3D2NBeamCrLinearTest
from NightlyTests import Simple3D2NBeamCrDynamicTest as T3D2NBeamCrDynamicTest
from NightlyTests import Simple2D2NBeamCrTest as T2D2NBeamCrTest
# Shell tests
### OLD Tests Start, will be removed soon, Philipp Bucher, 31.01.2018 |---
from NightlyTests import ShellQ4ThickBendingRollUpTests as TShellQ4ThickBendingRollUpTests
from NightlyTests import ShellQ4ThickDrillingRollUpTests as TShellQ4ThickDrillingRollUpTests
from NightlyTests import ShellQ4ThickOrthotropicLaminateLinearStaticTests as TShellQ4ThickOrthotropicLaminateLinearStaticTests
from NightlyTests import ShellT3ThinBendingRollUpTests as TShellT3ThinBendingRollUpTests
from NightlyTests import ShellT3ThinDrillingRollUpTests as TShellT3ThinDrillingRollUpTests
from NightlyTests import ShellT3IsotropicScordelisTests as TShellT3IsotropicScordelisTests
from NightlyTests import ShellT3ThinOrthotropicLaminateLinearStaticTests as TShellT3ThinOrthotropicLaminateLinearStaticTests
from NightlyTests import ShellT3ThickLinearStaticTests as TShellT3ThickLinearStaticTests
from NightlyTests import ShellT3ThickNonLinearStaticTests as TShellT3ThickNonLinearStaticTests
from NightlyTests import ShellT3ThickLinearDynamicTests as TShellT3ThickLinearDynamicTests
from NightlyTests import ShellT3ThickNonLinearDynamicTests as TShellT3ThickNonLinearDynamicTests
from NightlyTests import ShellT3ThickOrthotropicLaminateLinearStaticTests as TShellT3ThickOrthotropicLaminateLinearStaticTests
from NightlyTests import ShellQ4ThinLinearStaticTests as TShellQ4ThinLinearStaticTests
from NightlyTests import ShellQ4ThinNonLinearStaticTests as TShellQ4ThinNonLinearStaticTests
from NightlyTests import ShellQ4ThinLinearDynamicTests as TShellQ4ThinLinearDynamicTests
from NightlyTests import ShellQ4ThinNonLinearDynamicTests as TShellQ4ThinNonLinearDynamicTests
from NightlyTests import ShellQ4ThinOrthotropicLaminateLinearStaticTests as TShellQ4ThinOrthotropicLaminateLinearStaticTests
### ---| OLD Tests End
# Shell tests
from NightlyTests import ShellT3IsotropicLinearStaticStructScordelisLoRoofTests as TShellT3IsotropicLinearStaticStructScordelisLoRoofTests
from NightlyTests import ShellT3AndQ4LinearStaticStructScordelisLoRoofTests as TShellT3AndQ4LinearStaticStructScordelisLoRoofTests
from NightlyTests import ShellT3AndQ4LinearStaticStructPinchedCylinderTests as TShellT3AndQ4LinearStaticStructPinchedCylinderTests
from NightlyTests import ShellT3AndQ4LinearStaticStructPinchedHemisphereTests as TShellT3AndQ4LinearStaticStructPinchedHemisphereTests
from NightlyTests import ShellT3AndQ4LinearStaticStructClampedCylinderOrthotropicTests as TShellT3AndQ4LinearStaticStructClampedCylinderOrthotropicTests
from NightlyTests import ShellT3AndQ4NonLinearStaticStructHingedCylRoofSnapthroughTests as TShellT3AndQ4NonLinearStaticStructHingedCylRoofSnapthroughTests
from NightlyTests import ShellT3AndQ4NonLinearStaticStructHingedCylRoofSnapthroughOrthotropicTests as TShellT3AndQ4NonLinearStaticStructHingedCylRoofSnapthroughOrthotropicTests
from NightlyTests import ShellT3AndQ4NonLinearDynamicStructOscillatingPlateTests as TShellT3AndQ4NonLinearDynamicStructOscillatingPlateTests
from NightlyTests import ShellT3AndQ4NonLinearDynamicStructOscillatingPlateLumpedTests as TShellT3AndQ4NonLinearDynamicStructOscillatingPlateLumpedTests
# CL tests
from NightlyTests import IsotropicDamageSimoJuPSTest    as TIsotropicDamageSimoJuPSTest

##### VALIDATION TESTS #####
# SPRISM tests
from ValidationTests import SprismPanTests              as TSprismPanTests
# Pendulus Tests with Solid Elements
from ValidationTests import PendulusTLTest              as TPendulusTLTest
from ValidationTests import PendulusULTest              as TPendulusULTest
# Pendulus Tests with Shell Elements
from ValidationTests import ShellT3AndQ4NonLinearDynamicStructPendulusTests as TShellT3AndQ4NonLinearDynamicStructPendulusTests
from ValidationTests import ShellT3AndQ4NonLinearDynamicStructPendulusLumpedTests as TShellT3AndQ4NonLinearDynamicStructPendulusLumpedTests
from ValidationTests import ShellT3AndQ4NonLinearDynamicUnstructPendulusTests as TShellT3AndQ4NonLinearDynamicUnstructPendulusTests
from ValidationTests import ShellT3AndQ4NonLinearDynamicUnstructPendulusLumpedTests as TShellT3AndQ4NonLinearDynamicUnstructPendulusLumpedTests
# Shell tests
from ValidationTests import ShellT3AndQ4LinearStaticUnstructScordelisLoRoofTests as TShellT3AndQ4LinearStaticUnstructScordelisLoRoofTests
from ValidationTests import ShellT3AndQ4LinearStaticUnstructUnstructPinchedCylinderTests as TShellT3AndQ4LinearStaticUnstructUnstructPinchedCylinderTests
from ValidationTests import ShellT3AndQ4LinearStaticUnstructPinchedHemisphereTests as TShellT3AndQ4LinearStaticUnstructPinchedHemisphereTests
from ValidationTests import ShellT3AndQ4LinearStaticUnstructClampedCylinderOrthotropicTests as TShellT3AndQ4LinearStaticUnstructClampedCylinderOrthotropicTests
from ValidationTests import ShellT3AndQ4NonLinearStaticUnstructHingedCylRoofSnapthroughTests as TShellT3AndQ4NonLinearStaticUnstructHingedCylRoofSnapthroughTests
from ValidationTests import ShellT3AndQ4NonLinearStaticUnstructHingedCylRoofSnapthroughOrthotropicTests as TShellT3AndQ4NonLinearStaticUnstructHingedCylRoofSnapthroughOrthotropicTests
from ValidationTests import ShellT3AndQ4NonLinearDynamicUnstructOscillatingPlateTests as TShellT3AndQ4NonLinearDynamicUnstructOscillatingPlateTests
from ValidationTests import ShellT3AndQ4NonLinearDynamicUnstructOscillatingPlateLumpedTests as TShellT3AndQ4NonLinearDynamicUnstructOscillatingPlateLumpedTests


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
    # These tests are executed by the continuous integration tool, so they have to be very fast!
    # Execution time << 1 sec on a regular PC !!!
    # If the tests in the smallSuite take too long then merging to master will not be possible!
    smallSuite = suites['small'] # These tests are executed by the continuous integration tool
    nightSuite = suites['nightly'] # These tests are executed in the nightly build

    ### Adding the self-contained tests
    # Constitutive Law tests
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TTestConstitutiveLaw]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TTestMassCalculation]))
    # Solids
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TTestPatchTestSmallStrain]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TTestPatchTestLargeStrain]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TTestQuadraticElements]))
    # Shells
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TTestPatchTestShells]))
    nightSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TTestPatchTestShellsStressRec])) # TODO should be in smallSuite but is too slow
    nightSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TTestPatchTestShellsOrthotropic])) # TODO should be in smallSuite but is too slow
    # Trusses
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TTestTruss3D2N]))
    # Beams
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TTestCrBeam3D2N]))
    nightSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TTestCrBeam2D2N])) # TODO should be in smallSuite but is too slow
    # Membranes
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TTestPatchTestFormfinding]))
    # Loading Conditions
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TTestLoadingConditionsPoint]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TTestLoadingConditionsLine]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TTestLoadingConditionsSurface]))
    # Nodal Damping
    nightSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TNodalDampingTests])) # TODO should be in smallSuite but is too slow
    # Multipoint Constraint
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TTestMultipointConstraints]))

    ### Adding Small Tests
    # Dynamic basic tests (leave these in the smallSuite to have the Exection script tested)
    smallSuite.addTest(TDynamicBossakTests('test_execution'))
    smallSuite.addTest(TDynamicNewmarkTests('test_execution'))

    ### Adding Nightly Tests
    # Basic moving mesh test
    nightSuite.addTest(TSimpleMeshMovingTest('test_execution'))
    # Patch test Small Displacements
    nightSuite.addTest(TSDTwoDShearQuaPatchTest('test_execution'))
    nightSuite.addTest(TSDTwoDShearTriPatchTest('test_execution'))
    nightSuite.addTest(TSDTwoDTensionQuaPatchTest('test_execution'))
    nightSuite.addTest(TSDTwoDTensionTriPatchTest('test_execution'))
    nightSuite.addTest(TSDThreeDShearHexaPatchTest('test_execution'))
    nightSuite.addTest(TSDThreeDShearTetraPatchTest('test_execution'))
    nightSuite.addTest(TSDThreeDTensionHexaPatchTest('test_execution'))
    nightSuite.addTest(TSDThreeDTensionTetraPatchTest('test_execution'))
    # Patch test Total Lagrangian
    nightSuite.addTest(TTLTwoDShearQuaPatchTest('test_execution'))
    nightSuite.addTest(TTLTwoDShearTriPatchTest('test_execution'))
    nightSuite.addTest(TTLTwoDTensionQuaPatchTest('test_execution'))
    nightSuite.addTest(TTLTwoDTensionTriPatchTest('test_execution'))
    nightSuite.addTest(TTLThreeDShearHexaPatchTest('test_execution'))
    nightSuite.addTest(TTLThreeDShearTetraPatchTest('test_execution'))
    nightSuite.addTest(TTLThreeDTensionHexaPatchTest('test_execution'))
    nightSuite.addTest(TTLThreeDTensionTetraPatchTest('test_execution'))
    # Patch test Updated Lagrangian
    nightSuite.addTest(TULTwoDShearQuaPatchTest('test_execution'))
    nightSuite.addTest(TULTwoDShearTriPatchTest('test_execution'))
    nightSuite.addTest(TULTwoDTensionQuaPatchTest('test_execution'))
    nightSuite.addTest(TULTwoDTensionTriPatchTest('test_execution'))
    nightSuite.addTest(TULThreeDShearHexaPatchTest('test_execution'))
    nightSuite.addTest(TULThreeDShearTetraPatchTest('test_execution'))
    nightSuite.addTest(TULThreeDTensionHexaPatchTest('test_execution'))
    nightSuite.addTest(TULThreeDTensionTetraPatchTest('test_execution'))
    # SPRISM tests
    nightSuite.addTest(TSprismMembranePatchTests('test_execution'))
    nightSuite.addTest(TSprismBendingPatchTests('test_execution'))
    # Membrane tests
    nightSuite.addTest(TFofi4PointTentnoCableTests('test_execution'))
    nightSuite.addTest(TFofi4PointTentCableTests('test_execution'))
    nightSuite.addTest(TMembraneQ4PointLoadTests('test_execution'))
    nightSuite.addTest(TMembraneQ4TrussPointLoadTests('test_execution'))
    # 2Node Element tests
    nightSuite.addTest(T3D2NTrussDynamicTest('test_execution'))
    nightSuite.addTest(T3D2NTrussLinearTest('test_execution'))
    nightSuite.addTest(T3D2NTrussTest('test_execution'))
    nightSuite.addTest(T3D2NBeamCrTest('test_execution'))
    nightSuite.addTest(T3D2NBeamCrLinearTest('test_execution'))
    nightSuite.addTest(T3D2NBeamCrDynamicTest('test_execution'))
    # Shell tests
    nightSuite.addTest(TShellT3IsotropicLinearStaticStructScordelisLoRoofTests('test_execution'))
    nightSuite.addTest(TShellT3AndQ4LinearStaticStructScordelisLoRoofTests('test_execution'))
    nightSuite.addTest(TShellT3AndQ4LinearStaticStructPinchedCylinderTests('test_execution'))
    nightSuite.addTest(TShellT3AndQ4LinearStaticStructPinchedHemisphereTests('test_execution'))
    # nightSuite.addTest(TShellT3AndQ4LinearStaticStructClampedCylinderOrthotropicTests('test_execution'))
    nightSuite.addTest(TShellT3AndQ4NonLinearStaticStructHingedCylRoofSnapthroughTests('test_execution'))
    # nightSuite.addTest(TShellT3AndQ4NonLinearStaticStructHingedCylRoofSnapthroughOrthotropicTests('test_execution'))
    # nightSuite.addTest(TShellT3AndQ4NonLinearDynamicStructOscillatingPlateTests('test_execution'))
    # nightSuite.addTest(TShellT3AndQ4NonLinearDynamicStructOscillatingPlateLumpedTests('test_execution'))
    # Constitutive Law tests
    # nightSuite.addTest(TIsotropicDamageSimoJuPSTest('test_execution')) # FIXME: Needs get up to date

    if (missing_external_dependencies == False):
        if (hasattr(KratosMultiphysics.ExternalSolversApplication, "FEASTSolver")):
            # Eigenvalues tests
            smallSuite.addTest(TEigenQ4Thick2x2PlateTests('test_execution'))
            smallSuite.addTest(TEigenTL3D8NCubeTests('test_execution'))
            nightSuite.addTest(TEigen3D3NThinCircleTests('test_execution'))
            # Harmonic analysis test
            smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([THarmonicAnalysisTests]))
            # Element damping test
            nightSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TSpringDamperElementTests])) # TODO should be in smallSuite but is too slow
        else:
            print("FEASTSolver solver is not included in the compilation of the External Solvers Application")

    nightSuite.addTests(smallSuite)

    ### Adding Validation Tests
    # For very long tests that should not be in nighly and you can use to validate
    validationSuite = suites['validation']
    # SPRISM tests
    # validationSuite.addTest(TSprismPanTests('test_execution')) # FIXME: Needs get up to date
    validationSuite.addTest(T2D2NBeamCrTest('test_execution')) # TODO should be in nightSuite but is too slow
    # Pendulus tests with Solid Elements
    validationSuite.addTest(TPendulusTLTest('test_execution'))
    validationSuite.addTest(TPendulusULTest('test_execution'))
    # Pendulus Tests with Shell Elements
    # validationSuite.addTest(TShellT3AndQ4NonLinearDynamicStructPendulusTests('test_execution'))
    # validationSuite.addTest(TShellT3AndQ4NonLinearDynamicStructPendulusLumpedTests('test_execution'))
    # validationSuite.addTest(TShellT3AndQ4NonLinearDynamicUnstructPendulusTests('test_execution'))
    # validationSuite.addTest(TShellT3AndQ4NonLinearDynamicUnstructPendulusLumpedTests('test_execution'))
    # Shell tests
    validationSuite.addTest(TShellT3AndQ4LinearStaticUnstructScordelisLoRoofTests('test_execution'))
    validationSuite.addTest(TShellT3AndQ4LinearStaticUnstructUnstructPinchedCylinderTests('test_execution'))
    validationSuite.addTest(TShellT3AndQ4LinearStaticUnstructPinchedHemisphereTests('test_execution'))
    # validationSuite.addTest(TShellT3AndQ4LinearStaticUnstructClampedCylinderOrthotropicTests('test_execution'))
    validationSuite.addTest(TShellT3AndQ4NonLinearStaticUnstructHingedCylRoofSnapthroughTests('test_execution'))
    # validationSuite.addTest(TShellT3AndQ4NonLinearStaticUnstructHingedCylRoofSnapthroughOrthotropicTests('test_execution'))
    # validationSuite.addTest(TShellT3AndQ4NonLinearDynamicUnstructOscillatingPlateTests('test_execution'))
    # validationSuite.addTest(TShellT3AndQ4NonLinearDynamicUnstructOscillatingPlateLumpedTests('test_execution'))

    ### OLD Shell Tests Start, will be removed soon, Philipp Bucher, 31.01.2018 |---
    # They have been moved to validation temporarily until they will be removed
    validationSuite.addTest(TShellT3ThinDrillingRollUpTests('test_execution'))
    validationSuite.addTest(TShellT3IsotropicScordelisTests('test_execution'))
    validationSuite.addTest(TShellQ4ThickBendingRollUpTests('test_execution'))
    validationSuite.addTest(TShellQ4ThickDrillingRollUpTests('test_execution'))
    validationSuite.addTest(TShellQ4ThickOrthotropicLaminateLinearStaticTests('test_execution'))
    validationSuite.addTest(TShellT3ThinBendingRollUpTests('test_execution'))
    validationSuite.addTest(TShellT3ThinOrthotropicLaminateLinearStaticTests('test_execution'))
    validationSuite.addTest(TShellT3ThickLinearStaticTests('test_execution'))
    validationSuite.addTest(TShellT3ThickNonLinearStaticTests('test_execution'))
    validationSuite.addTest(TShellT3ThickLinearDynamicTests('test_execution'))
    validationSuite.addTest(TShellT3ThickNonLinearDynamicTests('test_execution'))
    validationSuite.addTest(TShellT3ThickOrthotropicLaminateLinearStaticTests('test_execution'))
    validationSuite.addTest(TShellQ4ThinLinearStaticTests('test_execution'))
    validationSuite.addTest(TShellQ4ThinNonLinearStaticTests('test_execution'))
    validationSuite.addTest(TShellQ4ThinLinearDynamicTests('test_execution'))
    validationSuite.addTest(TShellQ4ThinNonLinearDynamicTests('test_execution'))
    validationSuite.addTest(TShellQ4ThinOrthotropicLaminateLinearStaticTests('test_execution'))
    ### ---| OLD Shell Tests End

    # Create a test suit that contains all the tests:
    allSuite = suites['all']
    allSuite.addTests(nightSuite) # already contains the smallSuite
    allSuite.addTests(validationSuite)

    return suites


if __name__ == '__main__':
    KratosUnittest.runTests(AssambleTestSuites())
