import KratosMultiphysics

import KratosMultiphysics.KratosUnittest as KratosUnittest
import os

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

class TestTetrahedralMeshOrientationCheck(KratosUnittest.TestCase):

    def test_TetrahedralMeshOrientationCheck(self):
        KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
        current_model = KratosMultiphysics.Model()

        model_part = current_model.CreateModelPart("Main")
        model_part_io = KratosMultiphysics.ModelPartIO(GetFilePath("auxiliar_files_for_python_unittest/mdpa_files/cube_few_elements"))
        model_part_io.ReadModelPart(model_part)

        check_mesh = KratosMultiphysics.TetrahedralMeshOrientationCheck(model_part, True)
        check_mesh.Execute()

        self.assertGreater(model_part.ProcessInfo[KratosMultiphysics.FLAG_VARIABLE], 0.0)

    def test_QuadraticTetrahedralMeshOrientationCheck(self):
        KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
        current_model = KratosMultiphysics.Model()

        model_part = current_model.CreateModelPart("Main")
        model_part_io = KratosMultiphysics.ModelPartIO(GetFilePath("auxiliar_files_for_python_unittest/mdpa_files/quadratic_cube_few_elements"))
        model_part_io.ReadModelPart(model_part)

        check_mesh = KratosMultiphysics.TetrahedralMeshOrientationCheck(model_part, True)
        check_mesh.Execute()

        self.assertGreater(model_part.ProcessInfo[KratosMultiphysics.FLAG_VARIABLE], 0.0)

if __name__ == '__main__':
    KratosUnittest.main()

