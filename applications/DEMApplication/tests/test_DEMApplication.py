import KratosMultiphysics.KratosUnittest as KratosUnittest

#import test_guis
import test_kinematic_constraints
import test_particle_creator_destructor
import test_wall_creator_destructor
import test_analytics
import test_glued_particles
import test_restart
import test_DEM_2D
import test_DEM_2D_contact
import test_DEM_3D_contact
import test_DEM_2D_constitutive_laws
import test_DEM_3D_constitutive_laws
import test_DEM_2D_restitution
import test_DEM_3D_restitution
import test_DEM_2D_continuum_vs_discontinuum
import test_DEM_3D_continuum_vs_discontinuum
import test_DEM_3D_continuum
import test_DEM_2D_inlet
import test_DEM_3D_inlet
import test_inlet
import test_DEM_2D_control_module
import test_post_process
import test_friction_decay
import test_forces_and_moments
import test_history_dependent_CLs
import test_clusters
import test_DEM_schemes
import test_random_variable
import test_DEM_search_tolerance
import test_DEM_search_flags
import test_erase_particles
import test_search_nodes
import test_dem_3d_parallel_bond_model
import test_dem_3d_smooth_joint_model
import test_moving_periodic_boundary
import test_servo_control
import test_properties_measure_utility
import sys
sys.path.append('DEM3D_chung_ooi_tests/test1_data')
sys.path.append('DEM3D_chung_ooi_tests/test2_data')
sys.path.append('DEM3D_chung_ooi_tests/test3_data')
sys.path.append('DEM3D_chung_ooi_tests/test4_data')
import Chung_Ooi_test_1
import Chung_Ooi_test_2
import Chung_Ooi_test_3
import Chung_Ooi_test_4

def AssembleTestSuites():

    ''' Populates the test suites to run. At least, it should populate the suites:
    "small", "nightly" and "all"

    Return
    ------

    suites: A dictionary of suites
        The set of suites with its test_cases added.
    '''
    suites = KratosUnittest.KratosSuites

    # Create a test suit with the selected tests (Small tests):
    smallSuite = suites['small']
    #smallSuite.addTest(test_guis.TestGUIs("test_GUIs_1"))
    #smallSuite.addTest(test_guis.TestGUIs("test_GUIs_2"))
    smallSuite.addTest(test_particle_creator_destructor.TestParticleCreatorDestructor("test_CreateSphericParticle1"))
    smallSuite.addTest(test_wall_creator_destructor.TestWallCreatorDestructor("test_CreateWallTriangle"))
    smallSuite.addTest(test_analytics.TestAnalytics("test_Analytics_1"))
    #smallSuite.addTest(test_analytics.TestAnalytics("test_Analytics_2"))
    #smallSuite.addTest(test_analytics.TestAnalytics("test_Analytics_3"))
    smallSuite.addTest(test_glued_particles.TestGluedParticles("test_Glued_Particles_1"))
    smallSuite.addTest(test_DEM_2D.TestDEM2D("test_DEM2D_1"))
    smallSuite.addTest(test_DEM_3D_contact.TestDEM3DContact("test_DEM3D_contact"))
    smallSuite.addTest(test_DEM_2D_constitutive_laws.DEM2DConstitutiveLaws("test_DEM2D_ConstitutiveLaws1"))
    smallSuite.addTest(test_DEM_2D_constitutive_laws.DEM2DConstitutiveLaws("test_DEM2D_ConstitutiveLaws2"))
    smallSuite.addTest(test_DEM_2D_constitutive_laws.DEM2DConstitutiveLaws("test_DEM2D_ConstitutiveLaws3"))
    smallSuite.addTest(test_DEM_3D_constitutive_laws.DEM3DConstitutiveLaws("test_DEM3D_ConstitutiveLaws1"))
    smallSuite.addTest(test_DEM_3D_constitutive_laws.DEM3DConstitutiveLaws("test_DEM3D_ConstitutiveLaws2"))
    smallSuite.addTest(test_DEM_3D_constitutive_laws.DEM3DConstitutiveLaws("test_DEM3D_ConstitutiveLaws3"))
    smallSuite.addTest(test_DEM_2D_inlet.TestDEM2DInlet("test_DEM2D_inlet"))
    smallSuite.addTest(test_DEM_3D_inlet.TestDEM3DInlet("test_DEM3D_inlet"))
    smallSuite.addTest(test_inlet.TestPieceWiseLinearDEMInlet("test_piecewise_linear_inlet"))
    smallSuite.addTest(test_DEM_2D_restitution.TestDEM2DRestitution("test_DEM2D_restitution_1"))
    smallSuite.addTest(test_DEM_2D_restitution.TestDEM2DRestitution("test_DEM2D_restitution_2"))
    smallSuite.addTest(test_DEM_3D_restitution.TestDEM3DRestitution("test_DEM3D_restitution_1"))
    smallSuite.addTest(test_DEM_3D_restitution.TestDEM3DRestitution("test_DEM3D_restitution_2"))
    smallSuite.addTest(test_DEM_3D_continuum.TestDEM3DContinuum("test_DEM3D_continuum"))
    smallSuite.addTest(test_DEM_2D_control_module.TestDEM2DControlModule("test_DEM2DControlModule"))
    smallSuite.addTest(test_post_process.TestPostProcess("test_gid_printing_many_results"))
    smallSuite.addTest(test_friction_decay.TestFrictionDecay("test_Friction_Decay"))
    smallSuite.addTest(test_forces_and_moments.TestExternalForcesAndMoments("test_ForcesAndMoments"))
    smallSuite.addTest(test_history_dependent_CLs.TestHistoryDependentCLs("test_HistoryDependentCLs"))
    smallSuite.addTest(test_clusters.TestClusters("test_clusters_1"))
    smallSuite.addTest(test_DEM_schemes.TestDEMSchemes("test_ForwardEuler"))
    smallSuite.addTest(test_DEM_schemes.TestDEMSchemes("test_Taylor"))
    smallSuite.addTest(test_DEM_schemes.TestDEMSchemes("test_Symplectic"))
    smallSuite.addTest(test_DEM_schemes.TestDEMSchemes("test_Verlet"))
    smallSuite.addTest(test_random_variable.TestRandomVariable("test_random_variable"))
    smallSuite.addTest(test_DEM_2D_continuum_vs_discontinuum.TestDEM2DContinuumVsDiscontinuum("test_DEM2D_continuum_vs_discontinuum"))
    smallSuite.addTest(test_DEM_3D_continuum_vs_discontinuum.TestDEM3DContinuumVsDiscontinuum("test_DEM3D_continuum_vs_discontinuum"))
    smallSuite.addTest(test_DEM_2D_contact.TestDEM2DContact("test_DEM2D_contact"))
    smallSuite.addTest(test_kinematic_constraints.TestKinematicConstraints("test_KinematicConstraints_1"))
    smallSuite.addTest(test_DEM_search_flags.TestDEM3DSearchFlag("test_DEM3D_search"))
    smallSuite.addTest(test_erase_particles.TestDEMEraseParticlesWithDelay("test_erase_particles_no_delay"))
    smallSuite.addTest(test_erase_particles.TestDEMEraseParticlesWithDelay("test_erase_particles_little_delay"))
    smallSuite.addTest(test_erase_particles.TestDEMEraseParticlesWithDelay("test_erase_particles_with_delay"))
    smallSuite.addTest(test_search_nodes.TestSearchNodes("test_SearchNodesInTargetModelPart"))
    smallSuite.addTest(test_dem_3d_parallel_bond_model.TestParallelBondModel("test_ParallelBondModel_1"))
    smallSuite.addTest(test_dem_3d_smooth_joint_model.TestSmoothJointModel("test_SmoothJointModel_1"))
    smallSuite.addTest(test_moving_periodic_boundary.TestMovingPeriodicBoundary("test_MovingPeriodicBoundary"))
    smallSuite.addTest(test_servo_control.TestServoControl("test_ServoControl"))
    smallSuite.addTest(test_properties_measure_utility.TestPropertiesMeasureUtility("test_PropertiesMeasureUtility"))

    # Create a test suit with the selected tests plus all small tests
    nightSuite = suites['nightly']
    nightSuite.addTest(test_restart.TestRestartOneBall("test_execution"))
    nightSuite.addTest(test_restart.TestRestartTwoBalls("test_execution"))
    nightSuite.addTest(test_DEM_search_tolerance.TestSearchTolerance("test_SearchA"))
    nightSuite.addTest(test_DEM_search_tolerance.TestSearchTolerance("test_SearchB"))
    nightSuite.addTest(test_DEM_search_tolerance.TestSearchTolerance("test_SearchC"))
    nightSuite.addTest(test_DEM_search_tolerance.TestSearchTolerance("test_SearchD"))
    nightSuite.addTest(Chung_Ooi_test_1.ChungOoiTest1("test_Run"))
    nightSuite.addTest(Chung_Ooi_test_2.ChungOoiTest2("test_Run"))
    nightSuite.addTest(Chung_Ooi_test_3.ChungOoiTest3("test_Run"))
    nightSuite.addTest(Chung_Ooi_test_4.ChungOoiTest4("test_Run"))

    nightSuite.addTests(smallSuite)

    # For very long tests that should not be in nightly and you can use to validate
    validationSuite = suites['validation']
    validationSuite.addTests(nightSuite)

    # Create a test suit that contains all the tests:
    allSuite = suites['all']
    allSuite.addTests(validationSuite)

    return suites


if __name__ == '__main__':
    KratosUnittest.runTests(AssembleTestSuites())
