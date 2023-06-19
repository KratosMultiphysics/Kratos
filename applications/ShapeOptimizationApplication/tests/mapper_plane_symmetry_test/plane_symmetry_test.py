# Import Kratos core and apps
import KratosMultiphysics as KM
import KratosMultiphysics.ShapeOptimizationApplication as KSO

# Additional imports
from KratosMultiphysics.KratosUnittest import TestCase
from KratosMultiphysics.ShapeOptimizationApplication.mapper_factory import CreateMapper


class PlaneSymmetryMapperTest(TestCase):

    def test_mapping(self):
        model = KM.Model()
        model_part = model.CreateModelPart("cylinder")
        model_part.ProcessInfo.SetValue(KM.DOMAIN_SIZE, 3)
        model_part.AddNodalSolutionStepVariable(KSO.DF1DX)
        model_part.AddNodalSolutionStepVariable(KSO.DF1DX_MAPPED)
        model_part.AddNodalSolutionStepVariable(KSO.SHAPE_UPDATE)

        model_part.CreateNewNode(1, 0,0,0)
        model_part.CreateNewNode(2, 1,0,0)
        model_part.CreateNewNode(3, 2,0,0)
        model_part.CreateNewNode(4, 3,0,0)
        model_part.CreateNewNode(5, 4,0,0)
        model_part.CreateNewNode(6, 5,0,0)
        model_part.CreateNewNode(7, 6,0,0)

        settings = KM.Parameters("""
        {
            "filter_radius"              : 3.0,
            "plane_symmetry"             : true,
            "plane_symmetry_settings"    : {
                "point" : [3.0, 0.0, 0.0],
                "normal": [1.0, 0.0, 0.0]
            }
        }""")
        mapper = CreateMapper(model_part, model_part, settings)

        model_part.Nodes[1].SetSolutionStepValue(KSO.DF1DX, [1.0, 0.5, 0.5])
        mapper.InverseMap(KSO.DF1DX, KSO.DF1DX_MAPPED)

        expected_dfdx_mapped = {
            1: [0.25,0.125,0.125],
            2: [0.16666666666666666,0.08333333333333333,0.08333333333333333],
            3: [0.08333333333333333,0.041666666666666664,0.041666666666666664],
            4: [0.0,0.0,0.0],
            5: [-0.08333333333333333,0.041666666666666664,0.041666666666666664],
            6: [-0.16666666666666666,0.08333333333333333,0.08333333333333333],
            7: [-0.25,0.125,0.125],
        }

        for node in model_part.Nodes:
            v = node.GetSolutionStepValue(KSO.DF1DX_MAPPED)
            self.assertVectorAlmostEqual(v, expected_dfdx_mapped[node.Id], 7)
            # print(f"{node.Id}: [{v[0]},{v[1]},{v[2]}],")

        mapper.Map(KSO.DF1DX_MAPPED, KSO.SHAPE_UPDATE)

        expected_update = {
            1: [0.19444444444444445,0.09722222222222222,0.09722222222222222],
            2: [0.14583333333333334,0.07291666666666667,0.07291666666666667],
            3: [0.08333333333333334,0.05092592592592593,0.05092592592592593],
            4: [0.0,0.037037037037037035,0.037037037037037035],
            5: [-0.08333333333333334,0.05092592592592593,0.05092592592592593],
            6: [-0.14583333333333334,0.07291666666666667,0.07291666666666667],
            7: [-0.19444444444444445,0.09722222222222222,0.09722222222222222],
        }

        for node in model_part.Nodes:
            v = node.GetSolutionStepValue(KSO.SHAPE_UPDATE)
            self.assertVectorAlmostEqual(v, expected_update[node.Id], 7)
            # print(f"{node.Id}: [{v[0]},{v[1]},{v[2]}],")

if __name__ == '__main__':
    KM.KratosUnittest.main()

