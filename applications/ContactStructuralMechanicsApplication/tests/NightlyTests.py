from contact_structural_mechanics_test_factory import ContactStructuralMechanicsTestFactory as TestFactory

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

class ALMTaylorPatchFrictionalTestContact(TestFactory):
    file_name = "ALM_frictional_contact_test_2D/taylor_patch_test"

class ALMPureFrictionalTestContact(TestFactory):
    file_name = "ALM_frictional_contact_test_2D/pure_friction_test"
