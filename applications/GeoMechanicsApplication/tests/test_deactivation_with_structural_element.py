import os

import KratosMultiphysics.KratosUnittest as KratosUnittest
import test_helper
import glob

import KratosMultiphysics.GeoMechanicsApplication

from KratosMultiphysics.GeoMechanicsApplication.geomechanics_analysis import GeoMechanicsAnalysis
class KratosGeoMechanicsDeactivationWithStructuralTest(KratosUnittest.TestCase):
    """
    This class contains a test for a multistage computation with deactivation of structural elements in the first stage
    """

    def test_deactivation_with_structural_element(self):
        """

        """
        test_name    = 'deactivation_with_structural_element'
        project_path = test_helper.get_file_path(test_name)

        file_pattern = os.path.join(project_path, "ProjectParameters_stage*.json")
        stage_files = glob.glob(file_pattern)
        print(stage_files)

        cwd = os.getcwd()
        os.chdir(project_path)
        model = KratosMultiphysics.Model()

        for parameter_file in stage_files:
            with open(parameter_file,'r') as parameter_file:
                parameters = KratosMultiphysics.Parameters(parameter_file.read())
            simulation = GeoMechanicsAnalysis(model,parameters)
            simulation.Run()

        os.chdir(cwd)

if __name__ == '__main__':
    KratosUnittest.main()
