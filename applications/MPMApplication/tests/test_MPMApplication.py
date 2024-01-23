# import Kratos
import KratosMultiphysics
import run_cpp_unit_tests

# Import Kratos "wrapper" for unittests
import KratosMultiphysics.KratosUnittest as KratosUnittest


# Import from Test Factories (with general analysis flows)
from mpm_test_factory import AxisSymmetricCircularPlate2DTriTest as TAxisSymmetricCircularPlate2DTriTest

from mpm_test_factory import BeamCantileverStaticLinearElasticPointLoad2DTriTest as TBeamCantileverStaticLinearElasticPointLoad2DTriTest
from mpm_test_factory import BeamCantileverStaticLinearElasticParticlePointLoad2DTriTest as TBeamCantileverStaticLinearElasticParticlePointLoad2DTriTest
from mpm_test_factory import BeamCantileverStaticLinearElasticLineLoad2DQuadTest as TBeamCantileverStaticLinearElasticLineLoad2DQuadTest
from mpm_test_factory import BeamCantileverStaticLinearElasticSurfaceLoad3DHexaTest as TBeamCantileverStaticLinearElasticSurfaceLoad3DHexaTest
from mpm_test_factory import BeamCantileverStaticHyperelasticSelfWeightLoad2DQuadTest as TBeamCantileverStaticHyperelasticSelfWeightLoad2DQuadTest
from mpm_test_factory import BeamCantileverLinearStaticHyperelasticSelfWeightLoad2DQuadTest as TBeamCantileverLinearStaticHyperelasticSelfWeightLoad2DQuadTest
from mpm_test_factory import BeamCantileverDynamicConsistentMassTest as TBeamCantileverDynamicConsistentMassTest
from mpm_test_factory import BeamCantileverDynamicHyperelasticUPTest as TBeamCantileverDynamicHyperelasticUPTest

from mpm_test_factory import CooksMembraneCompressibleTest as TCooksMembraneCompressibleTest
from mpm_test_factory import CooksMembraneUPCompressibleTest as TCooksMembraneUPCompressibleTest
from mpm_test_factory import CooksMembraneUPIncompressibleTest as TCooksMembraneUPIncompressibleTest

from mpm_test_factory import CLLinearElastic3DQuadTest as TCLLinearElastic3DQuadTest
from mpm_test_factory import CLDispNewtonianFluidTest as TCLDispNewtonianFluidTest

from mpm_test_factory import GravityApplicationTest as TGravityApplicationTest
from mpm_test_factory import GravityTimeStepTableTest as TGravityTimeStepTableTest

from mpm_test_factory import PenaltyImpositionBeamCantileverStaticHyperelasticSelfWeightLoad2DQuadTest as TPenaltyImpositionBeamCantileverStaticHyperelasticSelfWeightLoad2DQuadTest

from mpm_test_factory import SlipBoundaryTest as TSlipBoundaryTest

from mpm_test_factory import ExplicitOscillatingPointUSLTest as TExplicitOscillatingPointUSLTest
from mpm_test_factory import ExplicitOscillatingPointUSFTest as TExplicitOscillatingPointUSFTest
from mpm_test_factory import ExplicitOscillatingPointMUSLTest as TExplicitOscillatingPointMUSLTest
from mpm_test_factory import ExplicitOscillatingPointCentralDifferenceTest as TExplicitOscillatingPointCentralDifferenceTest
from mpm_test_factory import ExplicitOscillatingPointYCompressibleTest as TExplicitOscillatingPointYCompressibleTest
from mpm_test_factory import ExplicitOscillatingPointGravityTest as TExplicitOscillatingPointGravityTest
from mpm_test_factory import ExplicitOscillatingPointTriTest as TExplicitOscillatingPointTriTest
from mpm_test_factory import ExplicitAxisymDiskTriCompressibleTest as TExplicitAxisymDiskTriCompressibleTest
from mpm_test_factory import ExplicitAxisymDiskQuadCompressibleTest as TExplicitAxisymDiskQuadCompressibleTest
from mpm_test_factory import Explicit3dHexCompressibleOscillatingPointTest as TExplicit3dHexCompressibleOscillatingPointTest
from mpm_test_factory import Explicit3dTetCompressibleOscillatingPointTest as TExplicit3dTetCompressibleOscillatingPointTest

from mpm_test_factory import PQMPMExplicitQuadTest as TPQMPMExplicitQuadTest
from mpm_test_factory import PQMPMExplicitTriTest as TPQMPMExplicitTriTest
from mpm_test_factory import PQMPMExplicitHexTest as TPQMPMExplicitHexTest

##### RESTART TESTS #####
from restart_tests import MPMRestartTestBeamStaticLineLoad2D  as TMPMRestartTestBeamStaticLineLoad2D
from restart_tests import MPMRestartTestDynamicCantilever2D    as TMPMRestartTestDynamicCantilever2D



# Import from Test Factories (with different analysis flows)
from test_generate_material_point_element   import TestGenerateMaterialPointElement   as TTestGenerateMPElement
from test_generate_material_point_condition import TestGenerateMaterialPointCondition as TTestGenerateMPCondition
from test_material_point_erase_process      import TestMaterialPointEraseProcess      as TTestMPEraseProcess
from test_search_material_point_element     import TestSearchMaterialPointElement     as TTestSearchMPElement
from test_search_material_point_condition   import TestSearchMaterialPointCondition   as TTestSearchMPCondition
from test_mpm_vtk_output_process            import TestMPMVtkOutputProcess            as TTestMPMVtkOutputProcess
from test_static_loading_conditions_point   import TestStaticLoadingConditionsPoint   as TTestStaticLoadingConditionsPoint
from test_static_loading_conditions_line    import TestStaticLoadingConditionsLine    as TTestStaticLoadingConditionsLine
from test_static_loading_conditions_surface import TestStaticLoadingConditionsSurface as TTestStaticLoadingConditionsSurface


def AssembleTestSuites():
    ''' Populates the test suites to run.

    Populates the test suites to run. At least, it should populate the suites:
    "small", "nightly" and "all"

    Return
    ------

    suites: A dictionary of suites
        The set of suites with its test_cases added.
    '''

    suites = KratosUnittest.KratosSuites

    ### Create a test suit with the selected tests (Small tests):
    ### These tests have to be very fast!
    ### Execution time << 1 sec on a regular PC !!!

    ## These tests are executed by the continuous integration tool
    smallSuite = suites['small']

    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TTestGenerateMPElement]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TTestGenerateMPCondition]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TTestMPEraseProcess]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TTestSearchMPElement]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TTestSearchMPCondition]))
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TTestMPMVtkOutputProcess]))

    # TODO: Look further into these three tests as they are still failing for AMatrix
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TTestStaticLoadingConditionsPoint]))    # FIXME:
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TTestStaticLoadingConditionsLine]))     # FIXME:
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([TTestStaticLoadingConditionsSurface]))  # FIXME:

    smallSuite.addTest(TCLLinearElastic3DQuadTest('test_execution'))
    smallSuite.addTest(TCLDispNewtonianFluidTest('test_execution'))
    smallSuite.addTest(TGravityApplicationTest('test_execution'))
    smallSuite.addTest(TGravityTimeStepTableTest('test_execution'))

    # TODO: Look further into this test as they are still failing for AMatrix
    smallSuite.addTest(TSlipBoundaryTest('test_execution')) # FIXME:

    ## These tests are executed in the nightly build
    nightSuite = suites['nightly']
    nightSuite.addTests(smallSuite)
    nightSuite.addTest(TAxisSymmetricCircularPlate2DTriTest('test_execution'))
    nightSuite.addTest(TBeamCantileverStaticLinearElasticPointLoad2DTriTest('test_execution'))
    nightSuite.addTest(TBeamCantileverStaticLinearElasticParticlePointLoad2DTriTest('test_execution'))
    nightSuite.addTest(TBeamCantileverStaticLinearElasticLineLoad2DQuadTest('test_execution'))
    nightSuite.addTest(TBeamCantileverStaticLinearElasticSurfaceLoad3DHexaTest('test_execution'))
    nightSuite.addTest(TBeamCantileverStaticHyperelasticSelfWeightLoad2DQuadTest('test_execution'))
    nightSuite.addTest(TCooksMembraneCompressibleTest('test_execution'))
    nightSuite.addTest(TCooksMembraneUPCompressibleTest('test_execution'))
    nightSuite.addTest(TCooksMembraneUPIncompressibleTest('test_execution'))
    nightSuite.addTest(TPenaltyImpositionBeamCantileverStaticHyperelasticSelfWeightLoad2DQuadTest('test_execution'))
    nightSuite.addTest(TBeamCantileverLinearStaticHyperelasticSelfWeightLoad2DQuadTest('test_execution'))
    nightSuite.addTest(TBeamCantileverDynamicConsistentMassTest('test_execution'))
    nightSuite.addTest(TBeamCantileverDynamicHyperelasticUPTest('test_execution'))


    nightSuite.addTest(TExplicitOscillatingPointUSLTest('test_execution'))
    nightSuite.addTest(TExplicitOscillatingPointUSFTest('test_execution'))
    nightSuite.addTest(TExplicitOscillatingPointMUSLTest('test_execution'))
    nightSuite.addTest(TExplicitOscillatingPointCentralDifferenceTest('test_execution'))
    nightSuite.addTest(TExplicitOscillatingPointYCompressibleTest('test_execution'))
    nightSuite.addTest(TExplicitOscillatingPointGravityTest('test_execution'))
    nightSuite.addTest(TExplicitOscillatingPointTriTest('test_execution'))
    nightSuite.addTest(TExplicitAxisymDiskTriCompressibleTest('test_execution'))
    nightSuite.addTest(TExplicitAxisymDiskQuadCompressibleTest('test_execution'))
    nightSuite.addTest(TExplicit3dHexCompressibleOscillatingPointTest('test_execution'))
    nightSuite.addTest(TExplicit3dTetCompressibleOscillatingPointTest('test_execution'))

    nightSuite.addTest(TPQMPMExplicitQuadTest('test_execution'))
    nightSuite.addTest(TPQMPMExplicitTriTest('test_execution'))
    nightSuite.addTest(TPQMPMExplicitHexTest('test_execution'))

    nightSuite.addTest(TMPMRestartTestBeamStaticLineLoad2D('test_execution'))
    nightSuite.addTest(TMPMRestartTestDynamicCantilever2D('test_execution'))

    ### Adding Validation Tests
    ## For very long tests that should not be in nighly and you can use to validate
    validationSuite = suites['validation']

    ### Create a test suit that contains all the tests:
    allSuite = suites['all']
    allSuite.addTests(nightSuite) # already contains smallSuite for visibility
    allSuite.addTests(validationSuite)

    return suites

if __name__ == '__main__':
    KratosMultiphysics.Logger.PrintInfo("Unittests", "\nRunning cpp unit tests ...")
    run_cpp_unit_tests.run()
    KratosMultiphysics.Logger.PrintInfo("Unittests", "Finished running cpp unit tests!")

    KratosMultiphysics.Logger.PrintInfo("Unittests", "\nRunning python tests ...")
    KratosUnittest.runTests(AssembleTestSuites())
    KratosMultiphysics.Logger.PrintInfo("Unittests", "Finished python tests!")
