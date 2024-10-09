# Kratos Imports
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics

# Standard Imports
import pathlib

def GetFullPathToFile(fileName):
    return pathlib.Path(__file__).absolute().parent / fileName

def GenerateModelPart(generate_moments=True):
    model = KratosMultiphysics.Model()
    model_part = model.CreateModelPart("Main")

    model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
    model_part.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)
    model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
    if generate_moments:
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ROTATION)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.MOMENT)

    model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VISCOSITY)
    KratosMultiphysics.ModelPartIO(
        str(GetFullPathToFile("test_files/mdpa_files/test_model_part_io_read"))
    ).ReadModelPart(model_part)

    model_part.GetNode(1).SetSolutionStepValue(KratosMultiphysics.REACTION_X, 10.0)
    model_part.GetNode(2).SetSolutionStepValue(KratosMultiphysics.REACTION_Z, 20.0)
    if (generate_moments):
        model_part.GetNode(3).SetSolutionStepValue(KratosMultiphysics.MOMENT_Z, 1.0)

    return model_part

class TestForceAndTorqueUtils(KratosUnittest.TestCase):
    def setUp(self):
        pass

    def test_force(self):
        model_part = GenerateModelPart()

        force = KratosMultiphysics.ForceAndTorqueUtils.SumForce(model_part, KratosMultiphysics.REACTION)

        self.assertAlmostEqual(force[0], 10.0)
        self.assertAlmostEqual(force[1], 0.0)
        self.assertAlmostEqual(force[2], 20.0)

    def test_with_moment(self):
        model_part = GenerateModelPart()

        Array3 = KratosMultiphysics.Array3
        reference_point = Array3([0.0, 0.0, 0.0])

        # Bare sum of forces and moments
        force, moment = KratosMultiphysics.ForceAndTorqueUtils.SumForceAndTorque(
            model_part,
            KratosMultiphysics.REACTION,
            KratosMultiphysics.MOMENT)

        self.assertAlmostEqual(force[0], 10.0)
        self.assertAlmostEqual(force[1], 0.0)
        self.assertAlmostEqual(force[2], 20.0)
        self.assertAlmostEqual(moment[0], 0.0)
        self.assertAlmostEqual(moment[1], 0.0)
        self.assertAlmostEqual(moment[2], 1.0)

        # Total reaction from the model part
        force, moment = KratosMultiphysics.ForceAndTorqueUtils.ComputeEquivalentForceAndTorque(
            model_part,
            reference_point,
            KratosMultiphysics.REACTION,
            KratosMultiphysics.MOMENT
        )

        self.assertAlmostEqual(force[0], 10.0)
        self.assertAlmostEqual(force[1], 0.0)
        self.assertAlmostEqual(force[2], 20.0)
        self.assertAlmostEqual(moment[0], 8.0)
        self.assertAlmostEqual(moment[1], -320.0)
        self.assertAlmostEqual(moment[2], 1.0)

    def test_without_moment(self):
        model_part = GenerateModelPart(generate_moments=False)

        Array3 = KratosMultiphysics.Array3
        reference_point = Array3([0.0, 0.0, 0.0])

        force, moment = KratosMultiphysics.ForceAndTorqueUtils.ComputeEquivalentForceAndTorque(
            model_part,
            reference_point,
            KratosMultiphysics.REACTION,
            KratosMultiphysics.MOMENT)

        self.assertAlmostEqual(force[0], 10.0)
        self.assertAlmostEqual(force[1], 0.0)
        self.assertAlmostEqual(force[2], 20.0)
        self.assertAlmostEqual(moment[0], 8.0)
        self.assertAlmostEqual(moment[1], -320.0)
        self.assertAlmostEqual(moment[2], 0.0)

if __name__ == "__main__":
    KratosUnittest.main()