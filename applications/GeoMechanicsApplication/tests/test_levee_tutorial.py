import os

import KratosMultiphysics as Kratos
import KratosMultiphysics.GeoMechanicsApplication as KratosGeo
import KratosMultiphysics.KratosUnittest as KratosUnittest
import test_helper

class KratosGeoMechanicsLeveeTutorialTests(KratosUnittest.TestCase):
    """
    This class contains a test case for the levee tutorial in the GeoMechanics application. 
    The results of the tutorial are verified against expected outcomes from reference simulations.
    The tutorial involves setting up a levee model, applying appropriate boundary conditions (including seepage), and running a simulation to observe the geomechanical behavior.
    """
    def test_levee_tutorial(self):
        # Define the path to the tutorial input file
        input_file = os.path.join(os.path.dirname(__file__), "levee_tutorial.json")

        # Run the tutorial using the test helper
        test_helper.run_tutorial(input_file)

        # Verify the results against expected outcomes
        test_helper.verify_tutorial_results(input_file)

        # Clean up any temporary files generated during the test
        test_helper.cleanup_tutorial_files(input_file)