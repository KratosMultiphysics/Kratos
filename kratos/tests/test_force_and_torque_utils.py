# Kratos Imports
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics

# Standard Imports
import pathlib

def GetFullPathToFile(fileName):
    return pathlib.Path(__file__).parent / fileName

class TestForceAndTorqueUtils(KratosUnittest.TestCase):
    def setUp(self):
        pass

    def GenerateModelPart(self, generate_moments=True):
        model = KratosMultiphysics.Model()
        model_part = model.CreateModelPart("Main")

        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)
        if generate_moments:
            model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ROTATION)
            model_part.AddNodalSolutionStepVariable(KratosMultiphysics.MOMENT)

        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VISCOSITY)
        KratosMultiphysics.ModelPartIO(
            str(GetFullPathToFile("auxiliar_files_for_python_unittest/mdpa_files/test_model_part_io_read"))
        ).ReadModelPart(model_part)

        model_part.GetNode(1).SetSolutionStepValue(KratosMultiphysics.REACTION_X, 10.0)
        model_part.GetNode(2).SetSolutionStepValue(KratosMultiphysics.REACTION_Z, 20.0)
        if (generate_moments):
            model_part.GetNode(3).SetSolutionStepValue(KratosMultiphysics.MOMENT_Z, 1.0)

        return model_part

    def test_with_moment(self):
        model_part = self.GenerateModelPart()

        Array3 = KratosMultiphysics.Array3
        referencePoint = Array3([0.0, 0.0, 0.0])
        force = Array3()
        moment = Array3()

        KratosMultiphysics.ForceAndTorqueUtils.SumForceAndTorque(
            model_part,
            KratosMultiphysics.Array3([0.0, 0.0, 0.0]),
            force,
            moment
        )

        self.assertAlmostEqual(force[0], 10.0)
        self.assertAlmostEqual(force[1], 0.0)
        self.assertAlmostEqual(force[2], 20.0)
        self.assertAlmostEqual(moment[0], 8.0)
        self.assertAlmostEqual(moment[1], -320.0)
        self.assertAlmostEqual(moment[2], 1.0)

    def test_without_moment(self):
        model_part = self.GenerateModelPart(generate_moments=False)

        Array3 = KratosMultiphysics.Array3
        referencePoint = Array3([0.0, 0.0, 0.0])
        force = Array3()
        moment = Array3()

        KratosMultiphysics.ForceAndTorqueUtils.SumForceAndTorque(
            model_part,
            KratosMultiphysics.Array3([0.0, 0.0, 0.0]),
            force,
            moment
        )

        self.assertAlmostEqual(force[0], 10.0)
        self.assertAlmostEqual(force[1], 0.0)
        self.assertAlmostEqual(force[2], 20.0)
        self.assertAlmostEqual(moment[0], 8.0)
        self.assertAlmostEqual(moment[1], -320.0)
        self.assertAlmostEqual(moment[2], 0.0)

if __name__ == "__main__":
    KratosUnittest.main()