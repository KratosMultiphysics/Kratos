from contact_structural_mechanics_test_factory import ContactStructuralMechanicsTestFactory as TestFactory

class ExplicitPenaltyThreeDSimplestPatchMatchingTestContact(TestFactory):
    file_name = "penalty_frictionless_contact_test_3D/explicit_3D_contact_simplest_patch_matching_test"

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

class ALMPureFrictionalTestContact(TestFactory):
    file_name = "ALM_frictional_contact_test_2D/pure_friction_test"

class ALMBasicFrictionTestContact(TestFactory):
    file_name = "ALM_frictional_contact_test_2D/basic_friction_test"

class ALMStaticEvolutionLoadFrictionTestContact(TestFactory):
    file_name = "ALM_frictional_contact_test_2D/static_evolution_load_test"

class ALMEvolutionLoadFrictionTestContact(TestFactory):
    file_name = "ALM_frictional_contact_test_2D/evolution_load_test"

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
