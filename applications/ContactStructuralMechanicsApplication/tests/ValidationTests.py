from contact_structural_mechanics_test_factory import ContactStructuralMechanicsTestFactory as TestFactory

class ALMTaylorPatchDynamicTestContact(TestFactory):
    file_name = "ALM_frictionless_contact_test_2D/taylor_patch_dynamic_test"

class ALMMeshMovingMatchingTestContact(TestFactory):
    file_name = "ALM_frictionless_contact_test_2D/mesh_moving_matching_test"

class ALMMeshMovingNotMatchingTestContact(TestFactory):
    file_name = "ALM_frictionless_contact_test_2D/mesh_moving_notmatching_test"

class ALMIroningTestContact(TestFactory):
    file_name = "ALM_frictionless_contact_test_2D/ironing_test"

class ALMIroningDieTestContact(TestFactory):
    file_name = "ALM_frictionless_contact_test_2D/ironing_die_test"

class LargeDisplacementPatchTestHexa(TestFactory):
    file_name = "mesh_tying_test/3D_contact_patch_test_large_disp_hexa"

class ALMLargeDisplacementPatchTestTetra(TestFactory):
    file_name = "ALM_frictionless_contact_test_3D/3D_contact_patch_test_large_disp_tetra"

class ALMLargeDisplacementPatchTestHexa(TestFactory):
    file_name = "ALM_frictionless_contact_test_3D/3D_contact_patch_test_large_disp_hexa"
    
class ALMMultiLayerContactTest(TestFactory):
    file_name = "ALM_frictionless_contact_test_3D/3D_multi_contact_test"

class ComponentsALMTaylorPatchDynamicTestContact(TestFactory):
    frictionless_by_components = True
    file_name = "ALM_frictionless_contact_test_2D/taylor_patch_dynamic_test"

class ComponentsALMMeshMovingMatchingTestContact(TestFactory):
    frictionless_by_components = True
    file_name = "ALM_frictionless_contact_test_2D/mesh_moving_matching_test"

class ComponentsALMMeshMovingNotMatchingTestContact(TestFactory):
    frictionless_by_components = True
    file_name = "ALM_frictionless_contact_test_2D/mesh_moving_notmatching_test"

class ComponentsALMLargeDisplacementPatchTestTetra(TestFactory):
    frictionless_by_components = True
    file_name = "ALM_frictionless_contact_test_3D/3D_contact_patch_test_large_disp_tetra"

class ComponentsALMLargeDisplacementPatchTestHexa(TestFactory):
    frictionless_by_components = True
    file_name = "ALM_frictionless_contact_test_3D/3D_contact_patch_test_large_disp_hexa"

class ComponentsALMMultiLayerContactTest(TestFactory):
    frictionless_by_components = True
    file_name = "ALM_frictionless_contact_test_3D/3D_multi_contact_test"
