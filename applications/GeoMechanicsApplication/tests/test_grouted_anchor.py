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


if __name__ == "__main__":
    KratosUnittest.main()
