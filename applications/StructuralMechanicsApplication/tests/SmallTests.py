import structural_mechanics_test_factory

# Do NOT remove these tests, they are there that the Execution script is being tested in the SmallTests
# (aka by the continuous integration tool before merging to master) !!!
class SimpleMeshMovingTest(structural_mechanics_test_factory.StructuralMechanicsTestFactory):
    file_name = "mesh_moving_test/simple_mesh_moving_test"
