import os
import sys

# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as KratosUtils
from KratosMultiphysics.testing.utilities import ReadModelPart

structural_mechanics_is_available = KratosUtils.CheckIfApplicationsAvailable("StructuralMechanicsApplication")
if structural_mechanics_is_available:
    import KratosMultiphysics.StructuralMechanicsApplication


def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)


class TestModelPartIOMPI(KratosUnittest.TestCase):
    def test_model_part_io_read_entity_data(self):
        # testing if the assignment of entity data works correctly in serial and MPI
        current_model = KratosMultiphysics.Model()
        model_part = current_model.CreateModelPart("Main")
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PARTITION_INDEX)
        model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] = 3
        ReadModelPart(GetFilePath("auxiliar_files_for_python_unittest/mdpa_files/cube_dummy_skin"), model_part)


if __name__ == '__main__':
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
    KratosUnittest.main()
