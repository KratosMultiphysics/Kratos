import os
from pathlib import Path
import KratosMultiphysics as Kratos
import KratosMultiphysics.GeoMechanicsApplication as KratosGeo
import KratosMultiphysics.KratosUnittest as KratosUnittest
import test_helper
import helper_utilities

class KratosGeoMechanicsLeveeTutorialTests(KratosUnittest.TestCase):
    """
    This class contains a test case for the levee tutorial in the GeoMechanics application. 
    The results of the tutorial are verified against expected outcomes from reference simulations.
    The tutorial involves setting up a levee model, applying appropriate boundary conditions (including seepage), and running a simulation to observe the geomechanical behavior.
    """
    def test_levee_seepage(self):
        # Define the path to the input files
        file_path = test_helper.get_file_path('test_levee_seepage')
        with KratosUnittest.WorkFolderScope(file_path, __file__):
            with open("ProjectParameters.json", "r") as f:
                project_parameters = Kratos.Parameters(f.read())
            helper_utilities.run_orchestrator(project_parameters)

if __name__ == '__main__':
    KratosUnittest.main()