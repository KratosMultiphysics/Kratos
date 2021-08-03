import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest

import os

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

class TestCadJsonInput(KratosUnittest.TestCase):

    def test_cad_json_input_read(self):
        cad_model = KratosMultiphysics.Model()
        cad_model_part = cad_model.CreateModelPart("CadModelPart")

        KratosMultiphysics.CadJsonInput(GetFilePath("auxiliar_files_for_python_unittest/cad_json_files/single_square")).ReadModelPart(cad_model_part)

        self.assertEqual(cad_model_part.NumberOfGeometries(), 8)
        self.assertTrue(cad_model_part.HasGeometry(1))
        self.assertFalse(cad_model_part.HasGeometry(10))
        self.assertTrue(cad_model_part.HasGeometry(55))

        # Check trimming edge index
        self.assertEqual(cad_model_part.GetGeometry(1).GetGeometryPart(2).Id, 2)

        background_geometry_index = KratosMultiphysics.Geometry.BACKGROUND_GEOMETRY_INDEX
        
        # Check unique indices
        self.assertNotEqual(cad_model_part.GetGeometry(1).GetGeometryPart(background_geometry_index).Id,
                            cad_model_part.GetGeometry(15).GetGeometryPart(background_geometry_index).Id)

if __name__ == '__main__':
    KratosUnittest.main()
