import KratosMultiphysics.KratosUnittest as KratosUnittest
import test_helper

import helper_utilities
import os

import KratosMultiphysics as Kratos


class KratosGeoMechanicsInterfaceTests(KratosUnittest.TestCase):
    def test_run(self):
        file_path = test_helper.get_file_path("grouted_anchor")

        with KratosUnittest.WorkFolderScope(
            rel_path_work_folder=file_path, file_path=__file__
        ):
            with open("ProjectParameters.json") as parameter_file:
                project_parameters = Kratos.Parameters(parameter_file.read())

            helper_utilities.run_orchestrator(project_parameters)
            run_stage_1_path = os.path.join(file_path, "Sheetpile_with_anchor.post.res")
            displacement_result = test_helper.get_nodal_variable_from_ascii(
                run_stage_1_path, "DISPLACEMENT"
            )

        characteristic_node_ids = [1, 2, 3, 4, 5]
        regression_displacement_results = [
            [-0.00076362, -0.104623, 0],
            [-2.55326e-05, -0.104993, 0],
            [0.000950976, -0.105002, 0],
            [-0.00256835, -0.108894, 0],
            [-0.000377458, -0.102321, 0],
        ]
        for node_id, reg_disp in zip(
            characteristic_node_ids, regression_displacement_results
        ):
            # result at time 1
            self.assertVectorAlmostEqual(displacement_result[1][node_id], reg_disp)


if __name__ == "__main__":
    KratosUnittest.main()
