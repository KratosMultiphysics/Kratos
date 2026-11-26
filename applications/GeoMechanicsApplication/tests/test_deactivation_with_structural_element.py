import os

import KratosMultiphysics.KratosUnittest as KratosUnittest
import test_helper

import KratosMultiphysics.GeoMechanicsApplication
from KratosMultiphysics.GeoMechanicsApplication.geomechanics_analysis import (
    GeoMechanicsAnalysis,
)


class KratosGeoMechanicsDeactivationWithStructuralTest(KratosUnittest.TestCase):
    """
    This class contains a test for a multistage computation with deactivation of structural elements in the first stage
    and the use of 'use_input_model_part' as input type for the second stage.
    """

    def test_deactivation_with_structural_element(self):
        test_name = "deactivation_with_structural_element"
        project_path = test_helper.get_file_path(test_name)

        stage_files = [os.path.join(project_path, filename) for filename in ["ProjectParameters_stage1.json", "ProjectParameters_stage2.json"]]

        cwd = os.getcwd()
        os.chdir(project_path)
        model = KratosMultiphysics.Model()

        # In the first stage, the truss is deactivated, so the expected y displacements are the same
        # for nodes 3 and 4. In the second stage, the truss is activated and since it's only connected to node 4,
        # the displacement of node 4 is half that of node 3 (the stiffness is doubled due to the truss).
        expected_y_displacements_node_3 = [0.05, -0.1]
        expected_y_displacements_node_4 = [0.05, -0.05]

        for (
            parameter_file_path,
            expected_y_displacement_node_3,
            expected_y_displacement_node_4,
        ) in zip(
            stage_files,
            expected_y_displacements_node_3,
            expected_y_displacements_node_4,
        ):
            with open(parameter_file_path, "r") as parameter_file:
                parameters = KratosMultiphysics.Parameters(parameter_file.read())
            simulation = GeoMechanicsAnalysis(model, parameters)
            simulation.Run()
            displacements = test_helper.get_displacement(simulation)
            self.assertAlmostEqual(
                displacements[2][1], expected_y_displacement_node_3, 4
            )
            self.assertAlmostEqual(
                displacements[3][1], expected_y_displacement_node_4, 4
            )

        os.chdir(cwd)


if __name__ == "__main__":
    KratosUnittest.main()
