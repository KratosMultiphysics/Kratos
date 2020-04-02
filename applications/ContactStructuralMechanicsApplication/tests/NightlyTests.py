from contact_structural_mechanics_test_factory import ContactStructuralMechanicsTestFactory as TestFactory

class SimplestPatchTestThreeDTriQuadMeshTying(TestFactory):
    file_name = "mesh_tying_test/3D_contact_simplest_patch_matching_triquad_test"

class SimplestPatchTestThreeDQuadTriMeshTying(TestFactory):
    file_name = "mesh_tying_test/3D_contact_simplest_patch_matching_quadtri_test"

class SimplePatchTestThreeDMeshTying(TestFactory):
    file_name = "mesh_tying_test/simple_patch_test_3D"

class ALMTwoDPatchComplexGeomTestContact(TestFactory):
    file_name = "ALM_frictionless_contact_test_2D/contact_patch_complex_geom_test"

class ALMTwoDPatchComplexGeomSlopeTestContact(TestFactory):
    file_name = "ALM_frictionless_contact_test_2D/contact_patch_complex_geom_slope_test"

class ALMSimplePatchTestContact(TestFactory):
    file_name = "ALM_frictionless_contact_test_2D/simple_patch_test"

class ALMSimpleSlopePatchTestContact(TestFactory):
    file_name = "ALM_frictionless_contact_test_2D/simple_slope_patch_test"

class ALMSimplePatchNotMatchingATestContact(TestFactory):
    file_name = "ALM_frictionless_contact_test_2D/simple_patch_notmatching_a_test"

class ALMSimplePatchNotMatchingBTestContact(TestFactory):
    file_name = "ALM_frictionless_contact_test_2D/simple_patch_notmatching_b_test"

class ALMThreeDSimplestPatchTestTriQuadContact(TestFactory):
    file_name = "ALM_frictionless_contact_test_3D/3D_contact_simplest_patch_matching_triquad_test"

class ALMThreeDSimplestPatchTestQuadTriContact(TestFactory):
    file_name = "ALM_frictionless_contact_test_3D/3D_contact_simplest_patch_matching_quadtri_test"

class ALMThreeDSimplestPatchMatchingAdaptativeTestContact(TestFactory):
    file_name = "ALM_frictionless_contact_test_3D/3D_contact_simplest_patch_matching_adaptative_test"

class ALMThreeDPatchComplexGeomTestContact(TestFactory):
    file_name = "ALM_frictionless_contact_test_3D/3D_contact_patch_complex_geom_test"

class ALMThreeDSimplestPatchMatchingSlopeTestContact(TestFactory):
    file_name = "ALM_frictionless_contact_test_3D/3D_contact_simplest_patch_matching_slope_test"

class ALMThreeDPatchMatchingTestContact(TestFactory):
    file_name = "ALM_frictionless_contact_test_3D/3D_contact_patch_matching_test"

class ALMThreeDPatchNotMatchingTestContact(TestFactory):
    file_name = "ALM_frictionless_contact_test_3D/3D_contact_patch_nonmatching_test"

class ALMTaylorPatchTestContact(TestFactory):
    file_name = "ALM_frictionless_contact_test_2D/taylor_patch_test"

class ALMHertzSimpleSphereTestContact(TestFactory):
    file_name = "ALM_frictionless_contact_test_2D/simple_hertz_sphere_plate_test"

class ALMHertzSimpleTestContact(TestFactory):
    file_name = "ALM_frictionless_contact_test_2D/hertz_simple_test"

class ALMHertzSphereTestContact(TestFactory):
    file_name = "ALM_frictionless_contact_test_2D/hertz_sphere_plate_test"

class ALMHertzCompleteTestContact(TestFactory):
    file_name = "ALM_frictionless_contact_test_2D/hertz_complete_test"

class ComponentsALMTwoDPatchComplexGeomTestContact(TestFactory):
    frictionless_by_components = True
    file_name = "ALM_frictionless_contact_test_2D/contact_patch_complex_geom_test"

class ComponentsALMTwoDPatchComplexGeomSlopeTestContact(TestFactory):
    frictionless_by_components = True
    file_name = "ALM_frictionless_contact_test_2D/contact_patch_complex_geom_slope_test"

class ComponentsALMSimplePatchTestContact(TestFactory):
    frictionless_by_components = True
    file_name = "ALM_frictionless_contact_test_2D/simple_patch_test"

class ComponentsALMSimpleSlopePatchTestContact(TestFactory):
    frictionless_by_components = True
    file_name = "ALM_frictionless_contact_test_2D/simple_slope_patch_test"

class ComponentsALMSimplePatchNotMatchingATestContact(TestFactory):
    frictionless_by_components = True
    file_name = "ALM_frictionless_contact_test_2D/simple_patch_notmatching_a_test"

class ComponentsALMSimplePatchNotMatchingBTestContact(TestFactory):
    frictionless_by_components = True
    file_name = "ALM_frictionless_contact_test_2D/simple_patch_notmatching_b_test"

class ComponentsALMThreeDSimplestPatchTestTriQuadContact(TestFactory):
    frictionless_by_components = True
    file_name = "ALM_frictionless_contact_test_3D/3D_contact_simplest_patch_matching_triquad_test"

class ComponentsALMThreeDSimplestPatchTestQuadTriContact(TestFactory):
    frictionless_by_components = True
    file_name = "ALM_frictionless_contact_test_3D/3D_contact_simplest_patch_matching_quadtri_test"

class ComponentsALMThreeDSimplestPatchMatchingAdaptativeTestContact(TestFactory):
    frictionless_by_components = True
    file_name = "ALM_frictionless_contact_test_3D/3D_contact_simplest_patch_matching_adaptative_test"

class ComponentsALMThreeDPatchComplexGeomTestContact(TestFactory):
    frictionless_by_components = True
    file_name = "ALM_frictionless_contact_test_3D/3D_contact_patch_complex_geom_test"

class ComponentsALMThreeDSimplestPatchMatchingSlopeTestContact(TestFactory):
    frictionless_by_components = True
    file_name = "ALM_frictionless_contact_test_3D/3D_contact_simplest_patch_matching_slope_test"

class ComponentsALMThreeDPatchMatchingTestContact(TestFactory):
    frictionless_by_components = True
    file_name = "ALM_frictionless_contact_test_3D/3D_contact_patch_matching_test"

class ComponentsALMThreeDPatchNotMatchingTestContact(TestFactory):
    frictionless_by_components = True
    file_name = "ALM_frictionless_contact_test_3D/3D_contact_patch_nonmatching_test"

class ComponentsALMTaylorPatchTestContact(TestFactory):
    frictionless_by_components = True
    file_name = "ALM_frictionless_contact_test_2D/taylor_patch_test"

class ComponentsALMHertzSimpleSphereTestContact(TestFactory):
    frictionless_by_components = True
    file_name = "ALM_frictionless_contact_test_2D/simple_hertz_sphere_plate_test"

class ComponentsALMHertzSimpleTestContact(TestFactory):
    frictionless_by_components = True
    file_name = "ALM_frictionless_contact_test_2D/hertz_simple_test"

class ComponentsALMHertzSphereTestContact(TestFactory):
    frictionless_by_components = True
    file_name = "ALM_frictionless_contact_test_2D/hertz_sphere_plate_test"

class ComponentsALMHertzCompleteTestContact(TestFactory):
    frictionless_by_components = True
    file_name = "ALM_frictionless_contact_test_2D/hertz_complete_test"

class ALMPureFrictionalTestContact(TestFactory):
    file_name = "ALM_frictional_contact_test_2D/pure_friction_test"

class ALMBasicFrictionTestContact(TestFactory):
    file_name = "ALM_frictional_contact_test_2D/basic_friction_test"

class ALMStaticEvolutionLoadFrictionTestContact(TestFactory):
    file_name = "ALM_frictional_contact_test_2D/static_evolution_load_test"

class ALMEvolutionLoadFrictionTestContact(TestFactory):
    file_name = "ALM_frictional_contact_test_2D/evolution_load_test"

class ThreeDSimplestPatchMatchingSlopeTestContact(TestFactory):
    file_name = "mpc_contact_tests/3D_contact_simplest_patch_matching_slope_test"

class ThreeDPatchMatchingTestContact(TestFactory):
    file_name = "mpc_contact_tests/3D_contact_patch_matching_test"

class ThreeDPatchNotMatchingTestContact(TestFactory):
    file_name = "mpc_contact_tests/3D_contact_patch_nonmatching_test"

class BeamAxilSimpleContactTest(TestFactory):
    file_name = "mpc_contact_tests/contact_beams_axil_hexa_simple_test"

class BeamAxilContactTest(TestFactory):
    file_name = "mpc_contact_tests/contact_beams_axil_hexa_test"

class BeamAxilTetraContactTest(TestFactory):
    file_name = "mpc_contact_tests/contact_beams_axil_test"

class BeamContactTest(TestFactory):
    file_name = "mpc_contact_tests/beam_contact_static_test"

class BeamContactWithFrictionTest(TestFactory):
    file_name = "mpc_contact_tests/beam_contact_static_with_friction_test"

class BeamContactWithTyingTest(TestFactory):
    file_name = "mpc_contact_tests/beam_contact_static_with_tying_test"

class PlateTest(TestFactory):
    file_name = "mpc_contact_tests/plate_test"
