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
            1: [0.16666666666666669,0.08333333333333334,0.08333333333333334],
            2: [0.11111111111111112,0.05555555555555556,0.05555555555555556],
            3: [0.05555555555555556,0.02777777777777778,0.02777777777777778],
            4: [0.0,0.0,0.0],
            5: [-0.05555555555555556,0.02777777777777778,0.02777777777777778],
            6: [-0.11111111111111112,0.05555555555555556,0.05555555555555556],
            7: [-0.16666666666666669,0.08333333333333334,0.08333333333333334],
        }

        for node in model_part.Nodes:
            v = node.GetSolutionStepValue(KSO.DF1DX_MAPPED)
            self.assertVectorAlmostEqual(v, expected_dfdx_mapped[node.Id], 7)
            # print(f"{node.Id}: [{v[0]},{v[1]},{v[2]}],")

        mapper.Map(KSO.DF1DX_MAPPED, KSO.SHAPE_UPDATE)

        expected_update = {
            1: [0.08641975308641976,0.04320987654320988,0.04320987654320988],
            2: [0.06481481481481483,0.03240740740740741,0.03240740740740741],
            3: [0.03703703703703704,0.02263374485596708,0.02263374485596708],
            4: [0.0,0.01646090534979424,0.01646090534979424],
            5: [-0.03703703703703704,0.02263374485596708,0.02263374485596708],
            6: [-0.06481481481481483,0.03240740740740741,0.03240740740740741],
            7: [-0.08641975308641975,0.043209876543209874,0.043209876543209874],
        }

        for node in model_part.Nodes:
            v = node.GetSolutionStepValue(KSO.SHAPE_UPDATE)
            self.assertVectorAlmostEqual(v, expected_update[node.Id], 7)
            # print(f"{node.Id}: [{v[0]},{v[1]},{v[2]}],")

if __name__ == '__main__':
    KM.KratosUnittest.main()

