# Import Kratos core and apps
import KratosMultiphysics as KM
import KratosMultiphysics.OptimizationApplication as KOA

# Additional imports
from KratosMultiphysics.KratosUnittest import TestCase

class SymmetryUtilitiesTest(TestCase):

    def test_mapping(self):
        with KM.KratosUnittest.WorkFolderScope(".", __file__):
            model = KM.Model()
            model_part = model.CreateModelPart("hexagon")
            model_part.ProcessInfo.SetValue(KM.DOMAIN_SIZE, 3)
            model_part.AddNodalSolutionStepVariable(KOA.CX)
            model_part.AddNodalSolutionStepVariable(KOA.CT)
            model_part_io = KM.ModelPartIO("hexagon")
            model_part_io.ReadModelPart(model_part)

        # plane symmetry test for scalar and vector field
        plane_symm_settings = KM.Parameters("""
        {
            "name" : "plane_x",
            "type" : "plane_symmetry",
            "settings":{
                "point" : [0.0,0.0,0.0],
                "normal": [1,0,0]
            }
        }
        """)
        plane_symm_util = KOA.SymmetryUtility("plane_x",model_part,plane_symm_settings)
        plane_symm_util.Initialize()

        model_part.Nodes[17].SetSolutionStepValue(KOA.CX, [1.0, 0.0, 1.0])
        model_part.Nodes[17].SetSolutionStepValue(KOA.CT, 1.0)
        plane_symm_util.ApplyOnVectorField(KOA.CX)
        plane_symm_util.ApplyOnScalarField(KOA.CT)
        v = model_part.Nodes[14].GetSolutionStepValue(KOA.CX)
        s = model_part.Nodes[14].GetSolutionStepValue(KOA.CT)
        self.assertVectorAlmostEqual(v, [-0.5, 0.0, 0.5], 7)
        self.assertAlmostEqual(s, 0.5, 7)
        # reset values
        model_part.Nodes[17].SetSolutionStepValue(KOA.CX, [1.0, 0.0, 1.0])
        model_part.Nodes[17].SetSolutionStepValue(KOA.CT, 1.0)
        model_part.Nodes[14].SetSolutionStepValue(KOA.CX, [0.0, 0.0, 0.0])
        model_part.Nodes[14].SetSolutionStepValue(KOA.CT, 0.0)
        
        # rotational symmetry test for scalar and vector field
        rot_symm_settings = KM.Parameters("""
        {
            "name" : "rot_z",
            "type" : "rotational_symmetry",
            "settings":{
                "axis": [0,0,1],
                "angle": 60,
                "point": [0,0,0]
            }
        }
        """)

        rot_symm_util = KOA.SymmetryUtility("rot_z",model_part,rot_symm_settings)
        rot_symm_util.Initialize()
        
        rot_symm_util.ApplyOnVectorField(KOA.CX)
        rot_symm_util.ApplyOnScalarField(KOA.CT)
        v = model_part.Nodes[13].GetSolutionStepValue(KOA.CX)
        s = model_part.Nodes[13].GetSolutionStepValue(KOA.CT)
        self.assertVectorAlmostEqual(v, [-0.0833333, 0.144338, 0.166667], 5)
        self.assertAlmostEqual(s, 0.166667, 5)



if __name__ == '__main__':
    KM.KratosUnittest.main()

