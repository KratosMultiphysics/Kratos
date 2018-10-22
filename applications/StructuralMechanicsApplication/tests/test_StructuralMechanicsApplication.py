# import Kratos
import KratosMultiphysics
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import run_cpp_unit_tests

# Import Kratos "wrapper" for unittests
import KratosMultiphysics.KratosUnittest as KratosUnittest

import subprocess

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
from test_constitutive_law import TestConstitutiveLaw as TTestConstitutiveLaw
# Processes test
from test_mass_calculation import TestMassCalculation as TTestMassCalculation
# Simple patch tests
from test_patch_test_small_strain import TestPatchTestSmallStrain as TTestPatchTestSmallStrain
from test_patch_test_small_strain_bbar import TestPatchTestSmallStrainBbar as TTestPatchTestSmallStrainBbar
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
# Nodal damping test
from test_nodal_damping import NodalDampingTests as TNodalDampingTests
# Spring damper element test
from test_spring_damper_element import SpringDamperElementTests as TSpringDamperElementTests
# Harmonic analysis tests
from test_harmonic_analysis import HarmonicAnalysisTests as THarmonicAnalysisTests
# Dynamic basic tests
from test_dynamic_schemes import FastDynamicSchemesTests as TFastDynamicSchemesTests
from test_dynamic_schemes import DynamicSchemesTests as TDynamicSchemesTests
# Eigenvalues Postprocessing Process test
from test_postprocess_eigenvalues_process import TestPostprocessEigenvaluesProcess as TTestPostprocessEigenvaluesProcess
# local-axis visualization tests
from test_local_axis_visualization import TestLocalAxisVisualization as TTestLocalAxisVisualization
# Test adjoint elements
from test_cr_beam_adjoint_element_3d2n import TestCrBeamAdjointElement as TTestCrBeamAdjointElement
from test_linear_thin_shell_adjoint_element_3d3n import TestShellThinAdjointElement3D3N as TTestShellThinAdjointElement3D3N
from test_truss_adjoint_element_3d2n import TestTrussAdjointElement as TTestTrussAdjointElement
from test_truss_adjoint_element_3d2n import TestTrussLinearAdjointElement as TTestTrussLinearAdjointElement
from test_adjoint_sensitity_analysis_beam_3d2n_structure import TestAdjointSensitivityAnalysisBeamStructure as TTestAdjointSensitivityAnalysisBeamStructure
from test_adjoint_sensitity_analysis_shell_3d3n_structure import TestAdjointSensitivityAnalysisShell3D3NStructure as TTestAdjointSensitivityAnalysisShell3D3NStructure
from test_adjoint_sensitity_analysis_truss_3d2n_structure import TestAdjointSensitivityAnalysisLinearTrussStructure as TTestAdjointSensitivityAnalysisLinearTrussStructure
from test_adjoint_sensitity_analysis_truss_3d2n_structure import TestAdjointSensitivityAnalysisNonLinearTrussStructure as TTestAdjointSensitivityAnalysisNonLinearTrussStructure

##### SMALL TESTS #####
# Basic moving mesh test (leave these in the smallSuite to have the Exection script tested)
from structural_mechanics_test_factory import SimpleMeshMovingTest as TSimpleMeshMovingTest

##### NIGHTLY TESTS #####
# Patch test Small Displacements
from structural_mechanics_test_factory import SDTwoDShearQuaPatchTest as TSDTwoDShearQuaPatchTest
from structural_mechanics_test_factory import SDTwoDShearTriPatchTest as TSDTwoDShearTriPatchTest
from structural_mechanics_test_factory import SDTwoDTensionQuaPatchTest as TSDTwoDTensionQuaPatchTest
from structural_mechanics_test_factory import SDTwoDTensionTriPatchTest as TSDTwoDTensionTriPatchTest
from structural_mechanics_test_factory import SDThreeDShearHexaPatchTest as TSDThreeDShearHexaPatchTest
from structural_mechanics_test_factory import SDThreeDShearTetraPatchTest as TSDThreeDShearTetraPatchTest
from structural_mechanics_test_factory import SDThreeDTensionHexaPatchTest as TSDThreeDTensionHexaPatchTest
from structural_mechanics_test_factory import SDThreeDTensionTetraPatchTest as TSDThreeDTensionTetraPatchTest
# Patch test Total Lagrangian
from structural_mechanics_test_factory import TLTwoDShearQuaPatchTest as TTLTwoDShearQuaPatchTest
from structural_mechanics_test_factory import TLTwoDShearTriPatchTest as TTLTwoDShearTriPatchTest
from structural_mechanics_test_factory import TLTwoDTensionQuaPatchTest as TTLTwoDTensionQuaPatchTest
from structural_mechanics_test_factory import TLTwoDTensionTriPatchTest as TTLTwoDTensionTriPatchTest
from structural_mechanics_test_factory import TLThreeDShearHexaPatchTest as TTLThreeDShearHexaPatchTest
from structural_mechanics_test_factory import TLThreeDShearTetraPatchTest as TTLThreeDShearTetraPatchTest
from structural_mechanics_test_factory import TLThreeDTensionHexaPatchTest as TTLThreeDTensionHexaPatchTest
from structural_mechanics_test_factory import TLThreeDTensionTetraPatchTest as TTLThreeDTensionTetraPatchTest
# Patch test Updated Lagrangian
from structural_mechanics_test_factory import ULTwoDShearQuaPatchTest as TULTwoDShearQuaPatchTest
from structural_mechanics_test_factory import ULTwoDShearTriPatchTest as TULTwoDShearTriPatchTest
from structural_mechanics_test_factory import ULTwoDTensionQuaPatchTest as TULTwoDTensionQuaPatchTest
from structural_mechanics_test_factory import ULTwoDTensionTriPatchTest as TULTwoDTensionTriPatchTest
from structural_mechanics_test_factory import ULThreeDShearHexaPatchTest as TULThreeDShearHexaPatchTest
from structural_mechanics_test_factory import ULThreeDShearTetraPatchTest as TULThreeDShearTetraPatchTest
from structural_mechanics_test_factory import ULThreeDTensionHexaPatchTest as TULThreeDTensionHexaPatchTest
from structural_mechanics_test_factory import ULThreeDTensionTetraPatchTest as TULThreeDTensionTetraPatchTest
# SPRISM tests
from structural_mechanics_test_factory import SprismMembranePatchTests as TSprismMembranePatchTests
from structural_mechanics_test_factory import SprismBendingPatchTests as TSprismBendingPatchTests
# Eigenvalues tests
from structural_mechanics_test_factory import EigenQ4Thick2x2PlateTests as TEigenQ4Thick2x2PlateTests
from structural_mechanics_test_factory import EigenTL3D8NCubeTests as TEigenTL3D8NCubeTests
from structural_mechanics_test_factory import Eigen3D3NThinCircleTests as TEigen3D3NThinCircleTests
# Membrane tests
from structural_mechanics_test_factory import Fofi4PointTentnoCableTests as TFofi4PointTentnoCableTests
from structural_mechanics_test_factory import Fofi4PointTentCableTests as TFofi4PointTentCableTests
from structural_mechanics_test_factory import MembraneQ4PointLoadTests as TMembraneQ4PointLoadTests
from structural_mechanics_test_factory import MembraneQ4TrussPointLoadTests as TMembraneQ4TrussPointLoadTests
# 2Node Element tests
from structural_mechanics_test_factory import Simple3D2NTrussTest as T3D2NTrussTest
from structural_mechanics_test_factory import Simple3D2NTrussLinearTest as T3D2NTrussLinearTest
from structural_mechanics_test_factory import Simple3D2NTrussDynamicTest as T3D2NTrussDynamicTest
from structural_mechanics_test_factory import Simple3D2NTrussLinearCompressionPlasticTest as T3D2NTrussLinearCompressionPlasticTest
from structural_mechanics_test_factory import Simple3D2NTrussLinearTensionPlasticTest as T3D2NTrussLinearTensionPlasticTest
from structural_mechanics_test_factory import Simple3D2NTrussNonLinearSnapthroughPlasticTest as T3D2NTrussNonLinearSnapthroughPlasticTest
from structural_mechanics_test_factory import Simple3D2NTrussNonLinearTensionPlasticTest as T3D2NTrussNonLinearTensionPlasticTest
from structural_mechanics_test_factory import Simple3D2NBeamCrTest as T3D2NBeamCrTest
from structural_mechanics_test_factory import Simple3D2NBeamCrLinearTest as T3D2NBeamCrLinearTest
from structural_mechanics_test_factory import Simple3D2NBeamCrDynamicTest as T3D2NBeamCrDynamicTest
from structural_mechanics_test_factory import Simple2D2NBeamCrTest as T2D2NBeamCrTest
# Shell tests
### OLD Tests Start, will be removed soon, Philipp Bucher, 31.01.2018 |---
from structural_mechanics_test_factory import ShellQ4ThickBendingRollUpTests as TShellQ4ThickBendingRollUpTests
from structural_mechanics_test_factory import ShellQ4ThickDrillingRollUpTests as TShellQ4ThickDrillingRollUpTests
from structural_mechanics_test_factory import ShellQ4ThickOrthotropicLaminateLinearStaticTests as TShellQ4ThickOrthotropicLaminateLinearStaticTests
from structural_mechanics_test_factory import ShellT3ThinBendingRollUpTests as TShellT3ThinBendingRollUpTests
from structural_mechanics_test_factory import ShellT3ThinDrillingRollUpTests as TShellT3ThinDrillingRollUpTests
from structural_mechanics_test_factory import ShellT3IsotropicScordelisTests as TShellT3IsotropicScordelisTests
from structural_mechanics_test_factory import ShellT3ThinOrthotropicLaminateLinearStaticTests as TShellT3ThinOrthotropicLaminateLinearStaticTests
from structural_mechanics_test_factory import ShellT3ThickLinearStaticTests as TShellT3ThickLinearStaticTests
from structural_mechanics_test_factory import ShellT3ThickNonLinearStaticTests as TShellT3ThickNonLinearStaticTests
from structural_mechanics_test_factory import ShellT3ThickLinearDynamicTests as TShellT3ThickLinearDynamicTests
from structural_mechanics_test_factory import ShellT3ThickNonLinearDynamicTests as TShellT3ThickNonLinearDynamicTests
from structural_mechanics_test_factory import ShellT3ThickOrthotropicLaminateLinearStaticTests as TShellT3ThickOrthotropicLaminateLinearStaticTests
from structural_mechanics_test_factory import ShellQ4ThinLinearStaticTests as TShellQ4ThinLinearStaticTests
from structural_mechanics_test_factory import ShellQ4ThinNonLinearStaticTests as TShellQ4ThinNonLinearStaticTests
from structural_mechanics_test_factory import ShellQ4ThinLinearDynamicTests as TShellQ4ThinLinearDynamicTests
from structural_mechanics_test_factory import ShellQ4ThinNonLinearDynamicTests as TShellQ4ThinNonLinearDynamicTests
from structural_mechanics_test_factory import ShellQ4ThinOrthotropicLaminateLinearStaticTests as TShellQ4ThinOrthotropicLaminateLinearStaticTests
### ---| OLD Tests End
# Shell tests
from structural_mechanics_test_factory import ShellT3IsotropicLinearStaticStructScordelisLoRoofTests as TShellT3IsotropicLinearStaticStructScordelisLoRoofTests
from structural_mechanics_test_factory import ShellT3AndQ4LinearStaticStructScordelisLoRoofTests as TShellT3AndQ4LinearStaticStructScordelisLoRoofTests
from structural_mechanics_test_factory import ShellT3AndQ4LinearStaticStructPinchedCylinderTests as TShellT3AndQ4LinearStaticStructPinchedCylinderTests
from structural_mechanics_test_factory import ShellT3AndQ4LinearStaticStructPinchedHemisphereTests as TShellT3AndQ4LinearStaticStructPinchedHemisphereTests
from structural_mechanics_test_factory import ShellT3AndQ4LinearStaticStructClampedCylinderOrthotropicTests as TShellT3AndQ4LinearStaticStructClampedCylinderOrthotropicTests
from structural_mechanics_test_factory import ShellT3AndQ4NonLinearStaticStructHingedCylRoofSnapthroughTests as TShellT3AndQ4NonLinearStaticStructHingedCylRoofSnapthroughTests
from structural_mechanics_test_factory import ShellT3AndQ4NonLinearStaticStructHingedCylRoofSnapthroughOrthotropicTests as TShellT3AndQ4NonLinearStaticStructHingedCylRoofSnapthroughOrthotropicTests
from structural_mechanics_test_factory import ShellT3AndQ4NonLinearDynamicStructOscillatingPlateTests as TShellT3AndQ4NonLinearDynamicStructOscillatingPlateTests
from structural_mechanics_test_factory import ShellT3AndQ4NonLinearDynamicStructOscillatingPlateLumpedTests as TShellT3AndQ4NonLinearDynamicStructOscillatingPlateLumpedTests
# CL tests
from structural_mechanics_test_factory import IsotropicDamageSimoJuPSTest    as TIsotropicDamageSimoJuPSTest
# Rigid test
from structural_mechanics_test_factory import RigidFaceTestWithImposeRigidMovementProcess as TRigidFaceTestWithImposeRigidMovementProcess

##### VALIDATION TESTS #####
# SPRISM tests
from structural_mechanics_test_factory import SprismPanTests              as TSprismPanTests
# Pendulus Tests with Solid Elements
from structural_mechanics_test_factory import PendulusTLTest              as TPendulusTLTest
from structural_mechanics_test_factory import PendulusULTest              as TPendulusULTest
# Pendulus Tests with Shell Elements
from structural_mechanics_test_factory import ShellT3AndQ4NonLinearDynamicStructPendulusTests as TShellT3AndQ4NonLinearDynamicStructPendulusTests
from structural_mechanics_test_factory import ShellT3AndQ4NonLinearDynamicStructPendulusLumpedTests as TShellT3AndQ4NonLinearDynamicStructPendulusLumpedTests
from structural_mechanics_test_factory import ShellT3AndQ4NonLinearDynamicUnstructPendulusTests as TShellT3AndQ4NonLinearDynamicUnstructPendulusTests
from structural_mechanics_test_factory import ShellT3AndQ4NonLinearDynamicUnstructPendulusLumpedTests as TShellT3AndQ4NonLinearDynamicUnstructPendulusLumpedTests
# Shell tests
from structural_mechanics_test_factory import ShellT3AndQ4LinearStaticUnstructScordelisLoRoofTests as TShellT3AndQ4LinearStaticUnstructScordelisLoRoofTests
from structural_mechanics_test_factory import ShellT3AndQ4LinearStaticUnstructUnstructPinchedCylinderTests as TShellT3AndQ4LinearStaticUnstructUnstructPinchedCylinderTests
from structural_mechanics_test_factory import ShellT3AndQ4LinearStaticUnstructPinchedHemisphereTests as TShellT3AndQ4LinearStaticUnstructPinchedHemisphereTests
from structural_mechanics_test_factory import ShellT3AndQ4LinearStaticUnstructClampedCylinderOrthotropicTests as TShellT3AndQ4LinearStaticUnstructClampedCylinderOrthotropicTests
from structural_mechanics_test_factory import ShellT3AndQ4NonLinearStaticUnstructHingedCylRoofSnapthroughTests as TShellT3AndQ4NonLinearStaticUnstructHingedCylRoofSnapthroughTests
from structural_mechanics_test_factory import ShellT3AndQ4NonLinearStaticUnstructHingedCylRoofSnapthroughOrthotropicTests as TShellT3AndQ4NonLinearStaticUnstructHingedCylRoofSnapthroughOrthotropicTests
from structural_mechanics_test_factory import ShellT3AndQ4NonLinearDynamicUnstructOscillatingPlateTests as TShellT3AndQ4NonLinearDynamicUnstructOscillatingPlateTests
from structural_mechanics_test_factory import ShellT3AndQ4NonLinearDynamicUnstructOscillatingPlateLumpedTests as TShellT3AndQ4NonLinearDynamicUnstructOscillatingPlateLumpedTests

##### RESTART TESTS #####
from restart_tests import TestSmallDisplacement2D4N  as TTestSmallDisplacement2D4N
from restart_tests import TestTotalLagrangian2D3N    as TTestTotalLagrangian2D3N
from restart_tests import TestUpdatedLagrangian3D8N  as TTestUpdatedLagrangian3D8N

##### RESPONSE_FUNCTION #####
from structural_response_function_test_factory import TestAdjointStrainEnergyResponseFunction as TTestAdjointStrainEnergyResponseFunction
from structural_response_function_test_factory import TestAdjointDisplacementResponseFunction as TTestAdjointDisplacementResponseFunction
from structural_response_function_test_factory import TestAdjointStressResponseFunction as TTestAdjointStressResponseFunction
from structural_response_function_test_factory import TestMassResponseFunction as TTestMassResponseFunction
from structural_response_function_test_factory import TestStrainEnergyResponseFunction as TTestStrainEnergyResponseFunction
from structural_response_function_test_factory import TestEigenfrequencyResponseFunction as TTestEigenfrequencyResponseFunction

def AssembleTestSuites():
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
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TTestPatchTestSmallStrainBbar]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TTestPatchTestLargeStrain]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TTestQuadraticElements]))
    # Shells
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TTestPatchTestShells]))
    nightSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TTestPatchTestShellsStressRec])) # TODO should be in smallSuite but is too slow
    nightSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TTestPatchTestShellsOrthotropic])) # TODO should be in smallSuite but is too slow
    # Membranes
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TTestPatchTestFormfinding]))
    # Trusses
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TTestTruss3D2N]))
    # Beams
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TTestCrBeam3D2N]))
    nightSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TTestCrBeam2D2N])) # TODO should be in smallSuite but is too slow
    # Loading Conditions
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TTestLoadingConditionsPoint]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TTestLoadingConditionsLine]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TTestLoadingConditionsSurface]))
    # Nodal Damping
    nightSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TNodalDampingTests])) # TODO should be in smallSuite but is too slow
    # Dynamic basic tests
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TFastDynamicSchemesTests]))
    # Eigenvalues Postprocessing Process test
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TTestPostprocessEigenvaluesProcess]))
    # local-axis visualization tests
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TTestLocalAxisVisualization]))
    # Adjoint Elements
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TTestCrBeamAdjointElement]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TTestShellThinAdjointElement3D3N]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TTestTrussAdjointElement]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TTestTrussLinearAdjointElement]))

    ### Adding Small Tests
    # Basic moving mesh test (leave these in the smallSuite to have the Exection script tested)
    smallSuite.addTest(TSimpleMeshMovingTest('test_execution'))
    # Basic restart test (leave these in the smallSuite to have the Exection script tested)
    smallSuite.addTest(TTestSmallDisplacement2D4N('test_execution'))
    smallSuite.addTest(TTestTotalLagrangian2D3N('test_execution'))
    smallSuite.addTest(TTestUpdatedLagrangian3D8N('test_execution'))

    ### Adding Nightly Tests
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
    nightSuite.addTest(T3D2NTrussLinearCompressionPlasticTest('test_execution'))
    nightSuite.addTest(T3D2NTrussLinearTensionPlasticTest('test_execution'))
    nightSuite.addTest(T3D2NTrussNonLinearSnapthroughPlasticTest('test_execution'))
    nightSuite.addTest(T3D2NTrussNonLinearTensionPlasticTest('test_execution'))
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
    nightSuite.addTest(TRigidFaceTestWithImposeRigidMovementProcess('test_execution'))

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

    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TTestAdjointSensitivityAnalysisBeamStructure]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TTestAdjointSensitivityAnalysisShell3D3NStructure]))
    nightSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TTestAdjointSensitivityAnalysisLinearTrussStructure]))
    nightSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TTestAdjointSensitivityAnalysisNonLinearTrussStructure]))

    nightSuite.addTest(TTestMassResponseFunction('test_execution'))
    nightSuite.addTest(TTestStrainEnergyResponseFunction('test_execution'))
    nightSuite.addTest(TTestEigenfrequencyResponseFunction('test_execution'))
    nightSuite.addTest(TTestAdjointStrainEnergyResponseFunction('test_execution'))
    nightSuite.addTest(TTestAdjointDisplacementResponseFunction('test_execution'))
    nightSuite.addTest(TTestAdjointStressResponseFunction('test_execution'))

    # Dynamic basic tests
    nightSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TDynamicSchemesTests]))

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
    # validationSuite.addTest(TShellQ4ThickDrillingRollUpTests('test_execution'))
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
    KratosMultiphysics.Logger.PrintInfo("Unittests", "\nRunning cpp unit tests ...")
    run_cpp_unit_tests.run()
    KratosMultiphysics.Logger.PrintInfo("Unittests", "Finished running cpp unit tests!")

    KratosMultiphysics.Logger.PrintInfo("Unittests", "\nRunning mpi python tests ...")
    try:
        import KratosMultiphysics.mpi as KratosMPI
        import KratosMultiphysics.MetisApplication as MetisApplication
        import KratosMultiphysics.TrilinosApplication as TrilinosApplication
        p = subprocess.Popen(["mpiexec", "-np", "2", "python3", "test_StructuralMechanicsApplication_mpi.py"], stdout=subprocess.PIPE)
        p.wait()
        KratosMultiphysics.Logger.PrintInfo("Unittests", "Finished mpi python tests!")
    except ImportError:
        KratosMultiphysics.Logger.PrintInfo("Unittests", "mpi is not available!")

    KratosMultiphysics.Logger.PrintInfo("Unittests", "\nRunning python tests ...")
    KratosUnittest.runTests(AssembleTestSuites())
    KratosMultiphysics.Logger.PrintInfo("Unittests", "Finished python tests!")
