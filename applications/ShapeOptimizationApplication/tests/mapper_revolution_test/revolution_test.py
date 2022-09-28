# Import Kratos core and apps
import KratosMultiphysics as KM
import KratosMultiphysics.ShapeOptimizationApplication as KSO

# Additional imports
from KratosMultiphysics.KratosUnittest import TestCase
from KratosMultiphysics.ShapeOptimizationApplication.mapper_factory import CreateMapper


class RevolutionMapperTest(TestCase):

    def test_mapping(self):
        with KM.KratosUnittest.WorkFolderScope(".", __file__):
            model = KM.Model()
            model_part = model.CreateModelPart("cylinder")
            model_part.ProcessInfo.SetValue(KM.DOMAIN_SIZE, 3)
            model_part.AddNodalSolutionStepVariable(KSO.DF1DX)
            model_part.AddNodalSolutionStepVariable(KSO.DF1DX_MAPPED)
            model_part.AddNodalSolutionStepVariable(KSO.SHAPE_UPDATE)
            model_part_io = KM.ModelPartIO("cylinder")
            model_part_io.ReadModelPart(model_part)

        settings = KM.Parameters("""
        {
            "filter_radius"              : 1.5,
            "revolution"             : true,
            "revolution_settings"    : {
                "point" : [0.0, 0.0, 0.0],
                "axis": [0.0, 0.0, 1.0]
            }
        }""")
        mapper = CreateMapper(model_part, model_part, settings)

        model_part.Nodes[1].SetSolutionStepValue(KSO.DF1DX, [1.0, 0.5, 0.5])
        mapper.InverseMap(KSO.DF1DX, KSO.DF1DX_MAPPED)

        expected_dfdx_mapped = {
            1: [0.12500000000136238,0.06250000000068119,0.06250000000068119],
            2: [0.041666666667120794,0.020833333333560397,0.020833333333560397],
            3: [0.008373412264039669,0.13950317547149493,0.06249999999931884],
            4: [-0.11662658773561976,0.07700317547166517,0.06249999999931884],
            5: [-0.12500000000136238,-0.0625000000006812,0.06250000000068119],
            6: [-0.008373412264039724,-0.13950317547149493,0.06249999999931884],
            7: [0.11662658773561974,-0.07700317547166521,0.06249999999931884],
            8: [0.0027911374214073977,0.046501058491511926,0.020833333333560397],
            9: [-0.03887552924605399,0.025667725157781214,0.020833333333560397],
            10: [-0.041666666667120794,-0.0208333333335604,0.020833333333560397],
            11: [-0.0027911374214074185,-0.04650105849151192,0.020833333333560397],
            12: [0.038875529246053975,-0.025667725157781235,0.020833333333560397],
        }

        for node in model_part.Nodes:
            v = node.GetSolutionStepValue(KSO.DF1DX_MAPPED)
            self.assertVectorAlmostEqual(v, expected_dfdx_mapped[node.Id], 7)
            # print(f"{node.Id}: [{v[0]},{v[1]},{v[2]}],")

        mapper.Map(KSO.DF1DX_MAPPED, KSO.SHAPE_UPDATE)

        expected_update = {
            1: [0.10416666666621256,0.052083333333106316,0.05208333333310628],
            2: [0.06250000000045412,0.03125000000022707,0.03125000000022707],
            3: [0.006977843553419639,0.11625264622713272,0.052083333333163076],
            4: [-0.09718882311375796,0.06416931289354391,0.052083333333163076],
            5: [-0.10416666666621256,-0.052083333333106316,0.05208333333310628],
            6: [-0.00697784355341964,-0.11625264622713273,0.052083333333163076],
            7: [0.09718882311375796,-0.06416931289354393,0.052083333333163076],
            8: [0.004186706132088288,0.06975158773688779,0.0312500000001703],
            9: [-0.058313293868763195,0.038501587736462034,0.0312500000001703],
            10: [-0.06250000000045412,-0.03125000000022707,0.03125000000022707],
            11: [-0.004186706132088288,-0.06975158773688779,0.0312500000001703],
            12: [0.058313293868763195,-0.038501587736462034,0.0312500000001703],
        }

        for node in model_part.Nodes:
            v = node.GetSolutionStepValue(KSO.SHAPE_UPDATE)
            self.assertVectorAlmostEqual(v, expected_update[node.Id], 7)
            # print(f"{node.Id}: [{v[0]},{v[1]},{v[2]}],")

if __name__ == '__main__':
    KM.KratosUnittest.main()

