from contact_structural_mechanics_test_factory import ContactStructuralMechanicsTestFactory as TestFactory

class SimplePatchTestTwoDMeshTying(TestFactory):
    file_name = "mesh_tying_test/simple_patch_test_2D"

class SimpleSlopePatchTestTwoDMeshTying(TestFactory):
    file_name = "mesh_tying_test/hyper_simple_slope_patch_test_2D"

class SimplestPatchTestThreeDMeshTying(TestFactory):
    file_name = "mesh_tying_test/3D_contact_simplest_patch_matching_test"

class SimplestPatchTestThreeDTriQuadMeshTying(TestFactory):
    file_name = "mesh_tying_test/3D_contact_simplest_patch_matching_triquad_test"

class SimplestPatchTestThreeDQuadTriMeshTying(TestFactory):
    file_name = "mesh_tying_test/3D_contact_simplest_patch_matching_quadtri_test"

class SimplePatchTestThreeDMeshTying(TestFactory):
    file_name = "mesh_tying_test/simple_patch_test_3D"
    
class ALMHyperSimplePatchTestContact(TestFactory):
    file_name = "ALM_frictionless_contact_test_2D/hyper_simple_patch_test"

class ALMHyperSimpleSlopePatchTestContact(TestFactory):
    file_name = "ALM_frictionless_contact_test_2D/hyper_simple_slope_patch_test"

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

class ALMThreeDSimplestPatchMatchingTestContact(TestFactory):
    file_name = "ALM_frictionless_contact_test_3D/3D_contact_simplest_patch_matching_test"

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

class ComponentsALMHyperSimpleTrianglePatchTestContact(TestFactory):
    frictionless_by_components = True
    file_name = "ALM_frictionless_contact_test_2D/hyper_simple_triangle_patch_test"

class ComponentsALMHyperSimplePatchTestContact(TestFactory):
    frictionless_by_components = True
    file_name = "ALM_frictionless_contact_test_2D/hyper_simple_patch_test"
    
class ComponentsALMHyperSimpleSlopePatchTestContact(TestFactory):
    frictionless_by_components = True
    file_name = "ALM_frictionless_contact_test_2D/hyper_simple_slope_patch_test"

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

class ComponentsALMThreeDSimplestPatchMatchingTestContact(TestFactory):
    frictionless_by_components = True
    file_name = "ALM_frictionless_contact_test_3D/3D_contact_simplest_patch_matching_test"

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

class ALMHyperSimplePatchFrictionalTestContact(TestFactory):
    file_name = "ALM_frictional_contact_test_2D/hyper_simple_patch_test"
