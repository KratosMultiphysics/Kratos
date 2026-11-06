import KratosMultiphysics.KratosUnittest as KratosUnittest
import test_helper

import helper_utilities
import os

import KratosMultiphysics as Kratos


class KratosGeoMechanicsInterfaceTests(KratosUnittest.TestCase):
    def test_run(self):
        file_path = test_helper.get_file_path("grouted_anchor")

        with open(
            os.path.join(file_path, "ProjectParameters.json"), "r"
        ) as parameter_file:
            project_parameters = Kratos.Parameters(parameter_file.read())

        cwd = os.getcwd()
        os.chdir(file_path)
        helper_utilities.run_orchestrator(project_parameters)
        os.chdir(cwd)

        run_stage_1_path = os.path.join(file_path, "Sheetpile_with_anchor.post.res")
        displacement_result = test_helper.get_nodal_variable_from_ascii(run_stage_1_path, "DISPLACEMENT")
        characteristic_node_ids = [1,2,3,4,5]
        regression_displacement_results=[[-0.00075417, -0.104971, 0],
                                        [-2.74607e-05, -0.105343, 0],
                                        [0.000989452, -0.105352, 0],
                                        [-0.00286768, -0.109109, 0],
                                        [7.98339e-05, -0.10327, 0]]
        for node_id, reg_disp in zip(characteristic_node_ids, regression_displacement_results):
            # result at time 1
            self.assertVectorAlmostEqual(displacement_result[1][node_id], reg_disp)


if __name__ == "__main__":
    KratosUnittest.main()
