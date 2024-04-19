
import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA

import KratosMultiphysics.KratosUnittest as kratos_unittest
from KratosMultiphysics.OptimizationApplication.responses import additive_manufacturing_responses

class TestOverHangResponseFunction(kratos_unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = Kratos.Model()
        cls.PrepareTest()
        cls.response_util.Initialize()

    @classmethod
    def PrepareTest(cls):
        cls.root_model_part = cls.model.CreateModelPart("root_mp")
        cls.evaluated_model_part = cls.root_model_part.CreateSubModelPart("evaluated_mp")
        cls.controlled_model_part = cls.root_model_part.CreateSubModelPart("controlled_mp")

        cls.response_settings = Kratos.Parameters("""{
                            "evaluated_objects": ["root_mp.evaluated_mp"],
                            "control_types": ["shape"],
                            "controlled_objects": ["root_mp"],
                            "print_direction": [0,0,1],
                            "max_angle": 45.0,
                            "heaviside_beta": 25.0,
                            "penalty_factor": 2.0,
                            "gradient_settings":{
                                "step_size" : 1e-6,
                                "gradient_mode": "finite_differencing"
                            }
                        }""")
        cls.response_util = additive_manufacturing_responses.MaxOverhangAngleResponseFunction("overhang",cls.response_settings,cls.model)

        cls.root_model_part.CreateNewNode(1, 0.0,0.5656854244,-0.8485281366)
        cls.root_model_part.CreateNewNode(2, 0.4,0.5656854244,-0.8485281366)
        cls.root_model_part.CreateNewNode(3, 0.0,0.8485281366,-0.5656854244)
        cls.root_model_part.CreateNewNode(4, 0.4,0.8485281366,-0.5656854244)
        cls.root_model_part.CreateNewNode(5, 0.0,0.2828427122,-0.5656854244)
        cls.root_model_part.CreateNewNode(6, 0.4,0.2828427122,-0.5656854244)
        cls.root_model_part.CreateNewNode(7, 0.0,0.5656854244,-0.2828427122)
        cls.root_model_part.CreateNewNode(8, 0.4,0.5656854244,-0.2828427122)


        properties = cls.root_model_part.CreateNewProperties(1)
        properties[Kratos.DENSITY] = 2.0
        cls.root_model_part.CreateNewElement("Element3D8N", 1, [1,2,4,3,5,6,8,7], properties)

        cls.root_model_part.CreateNewCondition("SurfaceCondition3D4N", 1, [1,3,4,2], properties)
        cls.root_model_part.CreateNewCondition("SurfaceCondition3D4N", 2, [1,5,7,3], properties)
        cls.root_model_part.CreateNewCondition("SurfaceCondition3D4N", 3, [3,7,8,4], properties)
        cls.root_model_part.CreateNewCondition("SurfaceCondition3D4N", 4, [2,4,8,6], properties)
        cls.root_model_part.CreateNewCondition("SurfaceCondition3D4N", 5, [1,2,6,5], properties)


        cls.evaluated_model_part.AddNodes([1,2,4,3,5,6,8,7])
        cls.evaluated_model_part.AddConditions([1,2,4,3,5])


        cls.controlled_model_part.AddNodes([1,2,4,3,5,6,8,7])
        cls.controlled_model_part.AddElements([1])


    def test_Value(self):
        self.assertAlmostEqual(self.response_util.CalculateValue(), 0.2, 5)

    def test_Sensitivities(self):
        self.response_util.CalculateGradientsForTypesAndObjects(["shape"],["root_mp"])

        node_3_sens = self.root_model_part.Nodes[3].GetSolutionStepValue(KratosOA.D_MAX_OVERHANG_ANGLE_D_X)
        node_4_sens = self.root_model_part.Nodes[4].GetSolutionStepValue(KratosOA.D_MAX_OVERHANG_ANGLE_D_X)
        node_5_sens = self.root_model_part.Nodes[5].GetSolutionStepValue(KratosOA.D_MAX_OVERHANG_ANGLE_D_X)
        node_6_sens = self.root_model_part.Nodes[6].GetSolutionStepValue(KratosOA.D_MAX_OVERHANG_ANGLE_D_X)


        self.assertVectorAlmostEqual(node_3_sens, [-0.1,1.9799,-1.83848], 5)
        self.assertVectorAlmostEqual(node_4_sens, [0.1,1.9799,-1.83848], 5)
        self.assertVectorAlmostEqual(node_5_sens, [-0.1,-1.9799,-1.83848], 5)
        self.assertVectorAlmostEqual(node_6_sens, [0.1,-1.9799,-1.83848], 5)


if __name__ == "__main__":
    kratos_unittest.main()