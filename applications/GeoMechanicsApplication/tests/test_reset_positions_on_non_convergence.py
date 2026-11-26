import os
import shutil

import KratosMultiphysics.KratosUnittest as KratosUnittest

# Import the Kratos applications that we need
import KratosMultiphysics.LinearSolversApplication
import KratosMultiphysics.StructuralMechanicsApplication
import KratosMultiphysics.GeoMechanicsApplication.run_geo_settlement as run_geo_settlement
from KratosMultiphysics.GeoMechanicsApplication.gid_output_file_reader import (
    GiDOutputFileReader,
)

import test_helper


class KratosGeoMechanicsResetPositionsOnNonConvergence(KratosUnittest.TestCase):
    """
    Two elements in series are stretched: a linear elastic element and a Mohr-Coulomb element with a
    tension cut-off. The second time step does not converge in its first cycle, so the state is reverted
    to the start of that step. Since the mesh is not moved, the nodal positions must not be changed by
    that revert (see issue #15007).
    """

    def setUp(self):
        super().setUp()

        self.test_root = test_helper.get_file_path("reset_positions_on_non_convergence")
        self.test_path = os.path.join(self.test_root, self.get_test_dir_name())

        shutil.rmtree(self.test_path, ignore_errors=True)
        os.makedirs(self.test_path)

        for filename in [
            "ProjectParameters.json",
            "MaterialParameters.json",
            "mesh.mdpa",
        ]:
            shutil.copy(
                os.path.join(self.test_root, filename),
                os.path.join(self.test_path, filename),
            )

    def get_test_dir_name(self):
        raise RuntimeError(
            "This base class does not provide a generic test directory name"
        )

    def check_results(self):
        output_data = GiDOutputFileReader().read_output_from(
            os.path.join(self.test_path, "reset_positions_on_non_convergence.post.res")
        )

        # The halved time step (t = 1.5) proves that the second step was cut back
        self.assertVectorAlmostEqual(
            GiDOutputFileReader.get_time_steps_from_first_valid_result(output_data),
            [1.0, 1.5, 2.0],
        )

        displacements = GiDOutputFileReader.nodal_values_at_time(
            "DISPLACEMENT", 2.0, output_data, [2, 3, 4, 5]
        )
        for displacement, expected_x in zip(displacements, [-0.2, -0.4, -0.4, -0.2]):
            self.assertAlmostEqual(displacement[0], expected_x, places=6)

        strains = GiDOutputFileReader.element_integration_point_values_at_time(
            "ENGINEERING_STRAIN_TENSOR", 2.0, output_data, [1, 2]
        )
        stresses = GiDOutputFileReader.element_integration_point_values_at_time(
            "CAUCHY_STRESS_TENSOR", 2.0, output_data, [1, 2]
        )
        for element_strains, element_stresses, expected_strain_xx in zip(
            strains, stresses, [-0.2, -0.2]
        ):
            for strain, stress in zip(element_strains, element_stresses):
                self.assertAlmostEqual(strain[0], expected_strain_xx, places=6)
                self.assertAlmostEqual(stress[0], -20.0, places=4)


class KratosGeoMechanicsResetPositionsOnNonConvergencePyRoute(
    KratosGeoMechanicsResetPositionsOnNonConvergence
):
    def get_test_dir_name(self):
        return "python"

    def test_non_converged_step_does_not_move_mesh(self):
        simulation = test_helper.run_kratos(self.test_path)

        self.check_results()
        for node in simulation._GetSolver().GetComputingModelPart().Nodes:
            self.assertVectorAlmostEqual(
                [node.X, node.Y, node.Z], [node.X0, node.Y0, node.Z0]
            )


class KratosGeoMechanicsResetPositionsOnNonConvergenceCppRoute(
    KratosGeoMechanicsResetPositionsOnNonConvergence
):
    def get_test_dir_name(self):
        return "cpp"

    def test_non_converged_step_does_not_move_mesh(self):
        status = run_geo_settlement.run_stages(
            self.test_path, ["ProjectParameters.json"]
        )
        self.assertEqual(status, 0)

        self.check_results()


if __name__ == "__main__":
    KratosUnittest.main()
