from __future__ import print_function, absolute_import, division

import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest

import KratosMultiphysics.kratos_utilities as kratos_utils

import os

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

class TestCadJsonInput(KratosUnittest.TestCase):

    def test_cad_json_input_read(self):
        cad_model = KratosMultiphysics.Model("CadModel")
        cad_model_part = cad_model.CreateModelPart("CadModelPart")
        cad_json_parameters = KratosMultiphysics.Parameters(GetFilePath("auxiliar_files_for_python_unittest/cad_json_files/single_square.json"))

        KratosMultiphysics.CadJsonInput(cad_json_parameters).ReadModelPart(cad_model_part)

if __name__ == '__main__':
    KratosUnittest.main()
