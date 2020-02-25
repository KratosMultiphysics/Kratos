from __future__ import print_function, absolute_import, division

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

        self.assertEqual(cad_model_part.NumberOfGeometries(), 5)
        self.assertTrue(cad_model_part.HasGeometry(1))
        self.assertFalse(cad_model_part.HasGeometry(10))

if __name__ == '__main__':
    KratosUnittest.main()
