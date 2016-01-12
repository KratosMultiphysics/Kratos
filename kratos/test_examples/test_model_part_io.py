import unittest
from KratosMultiphysics import *

class TestModelPartIO(unittest.TestCase):

    def test_model_part_io_read_model_part(self):
        model_part = ModelPart("Main")
        model_part.AddNodalSolutionStepVariable(DISPLACEMENT)
        model_part.AddNodalSolutionStepVariable(VISCOSITY)
        model_part_io = ModelPartIO("test_model_part_io")
        model_part_io.ReadModelPart(model_part)

        self.assertEqual(model_part.NumberOfSubModelParts(), 0)

        self.assertEqual(model_part.NumberOfNodes(), 6)
        self.assertEqual(model_part.NumberOfProperties(), 1)
        self.assertEqual(model_part.NumberOfElements(), 4)
        self.assertEqual(model_part.NumberOfConditions(), 5)

        self.assertTrue(model_part.GetNode(1).IsFixed(DISPLACEMENT_X))
        self.assertTrue(model_part.GetNode(2).IsFixed(DISPLACEMENT_X))
        self.assertFalse(model_part.GetNode(3).IsFixed(DISPLACEMENT_X))
        self.assertFalse(model_part.GetNode(972).IsFixed(DISPLACEMENT_X))
        self.assertTrue(model_part.GetNode(973).IsFixed(DISPLACEMENT_X))
        self.assertTrue(model_part.GetNode(974).IsFixed(DISPLACEMENT_X))

        self.assertTrue(model_part.GetNode(1).IsFixed(DISPLACEMENT_Y))
        self.assertTrue(model_part.GetNode(2).IsFixed(DISPLACEMENT_Y))
        self.assertFalse(model_part.GetNode(3).IsFixed(DISPLACEMENT_Y))
        self.assertFalse(model_part.GetNode(972).IsFixed(DISPLACEMENT_Y))
        self.assertTrue(model_part.GetNode(973).IsFixed(DISPLACEMENT_Y))
        self.assertTrue(model_part.GetNode(974).IsFixed(DISPLACEMENT_Y))

        self.assertTrue(model_part.GetNode(1).IsFixed(DISPLACEMENT_Z))
        self.assertTrue(model_part.GetNode(2).IsFixed(DISPLACEMENT_Z))
        self.assertFalse(model_part.GetNode(3).IsFixed(DISPLACEMENT_Z))
        self.assertFalse(model_part.GetNode(972).IsFixed(DISPLACEMENT_Z))
        self.assertTrue(model_part.GetNode(973).IsFixed(DISPLACEMENT_Z))
        self.assertTrue(model_part.GetNode(974).IsFixed(DISPLACEMENT_Z))

        self.assertEqual(model_part.GetNode(1).GetSolutionStepValue(DISPLACEMENT_X),0.1)
        self.assertEqual(model_part.GetNode(2).GetSolutionStepValue(DISPLACEMENT_X),0.2)
        self.assertEqual(model_part.GetNode(3).GetSolutionStepValue(DISPLACEMENT_X),0.0)
        self.assertEqual(model_part.GetNode(972).GetSolutionStepValue(DISPLACEMENT_X),0.0)
        self.assertEqual(model_part.GetNode(973).GetSolutionStepValue(DISPLACEMENT_X),0.0)
        self.assertEqual(model_part.GetNode(974).GetSolutionStepValue(DISPLACEMENT_X),0.0)

        self.assertEqual(model_part.GetNode(1).GetSolutionStepValue(DISPLACEMENT_Y),0.0)
        self.assertEqual(model_part.GetNode(2).GetSolutionStepValue(DISPLACEMENT_Y),0.0)
        self.assertEqual(model_part.GetNode(3).GetSolutionStepValue(DISPLACEMENT_Y),0.0)
        self.assertEqual(model_part.GetNode(972).GetSolutionStepValue(DISPLACEMENT_Y),0.0)
        self.assertEqual(model_part.GetNode(973).GetSolutionStepValue(DISPLACEMENT_Y),0.000973)
        self.assertEqual(model_part.GetNode(974).GetSolutionStepValue(DISPLACEMENT_Y),0.000974)

        self.assertEqual(model_part.GetNode(1).GetSolutionStepValue(DISPLACEMENT_Z),0.0)
        self.assertEqual(model_part.GetNode(2).GetSolutionStepValue(DISPLACEMENT_Z),0.0)
        self.assertEqual(model_part.GetNode(3).GetSolutionStepValue(DISPLACEMENT_Z),0.0)
        self.assertEqual(model_part.GetNode(972).GetSolutionStepValue(DISPLACEMENT_Z),0.0)
        self.assertEqual(model_part.GetNode(973).GetSolutionStepValue(DISPLACEMENT_Z),0.0)
        self.assertEqual(model_part.GetNode(974).GetSolutionStepValue(DISPLACEMENT_Z),0.0)

        self.assertEqual(model_part.GetNode(1).GetSolutionStepValue(VISCOSITY),0.01)
        self.assertEqual(model_part.GetNode(2).GetSolutionStepValue(VISCOSITY),0.01)
        self.assertEqual(model_part.GetNode(3).GetSolutionStepValue(VISCOSITY),0.0)
        self.assertEqual(model_part.GetNode(972).GetSolutionStepValue(VISCOSITY),0.0)
        self.assertEqual(model_part.GetNode(973).GetSolutionStepValue(VISCOSITY),0.01)
        self.assertEqual(model_part.GetNode(974).GetSolutionStepValue(VISCOSITY),0.01)

    #def test_model_part_io_properties_block(self):
    #    model_part = ModelPart("Main")
    #    model_part_io = ModelPartIO("test_model_part_io")
    #    model_part_io.ReadProperties(model_part.Properties)

if __name__ == '__main__':
    unittest.main()