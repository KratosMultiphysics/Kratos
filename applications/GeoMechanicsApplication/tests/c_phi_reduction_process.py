import os
import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.ConstitutiveLawsApplication as KratosConstitutiveLaws
import KratosMultiphysics.GeoMechanicsApplication.context_managers as context_managers
import KratosMultiphysics.GeoMechanicsApplication.geomechanics_analysis as analysis
from KratosMultiphysics.GeoMechanicsApplication.gid_output_file_reader import GiDOutputFileReader

import test_helper

class KratosGeoMechanicsCPhiReductionProcess(KratosUnittest.TestCase):
    """
    This class contains benchmark tests which are checked with the original solution
    """

    @staticmethod
    def _run_stages_with_common_parameters(case_path, common_path, stages):
        model = Kratos.Model()
        with context_managers.set_cwd_to(case_path):
            for stage_number, stage_index in enumerate(stages):
                case_parameters_file = os.path.join(
                    case_path, f"ProjectParameters_stage{stage_index}.json"
                )
                project_parameters_file = (
                    case_parameters_file
                    if os.path.exists(case_parameters_file)
                    else os.path.join(common_path, f"ProjectParameters_stage{stage_index}.json")
                )
                with open(project_parameters_file, "r") as parameter_file:
                    stage_parameters = Kratos.Parameters(parameter_file.read())
                model_import_settings = stage_parameters["solver_settings"]["model_import_settings"]
                if stage_number == 0 and model_import_settings["input_type"].GetString() == "use_input_model_part":
                    mesh_file = os.path.relpath(os.path.join(common_path, "mesh"), case_path)
                    model_import_settings["input_type"].SetString("mdpa")
                    model_import_settings.AddEmptyValue("input_filename").SetString(mesh_file)
                elif stage_number > 0:
                    model_import_settings["input_type"].SetString("use_input_model_part")

                simulation = analysis.GeoMechanicsAnalysis(model, stage_parameters)
                simulation.Run()

                if stage_index == 1 and os.path.exists(case_parameters_file):
                    KratosGeoMechanicsCPhiReductionProcess._assert_initial_state_is_elastic(simulation)

    @staticmethod
    def _assert_initial_state_is_elastic(simulation, tolerance=1.0e-10):
        plastic_strains = test_helper.get_on_integration_points(
            simulation, KratosConstitutiveLaws.PLASTIC_STRAIN_VECTOR
        )
        plastic_dissipations = test_helper.get_on_integration_points(
            simulation, Kratos.PLASTIC_DISSIPATION
        )

        max_plastic_strain = max(
            sum(component * component for component in value) ** 0.5
            for element_values in plastic_strains
            for value in element_values
        )
        max_plastic_dissipation = max(
            abs(value)
            for element_values in plastic_dissipations
            for value in element_values
        )

        if max_plastic_strain > tolerance or max_plastic_dissipation > tolerance:
            raise AssertionError(
                "The full-strength initial state is not elastic: "
                f"max ||plastic strain|| = {max_plastic_strain}, "
                f"max plastic dissipation = {max_plastic_dissipation}."
            )

    @staticmethod
    def _get_displacement_vectors(case_path):
        common_path = os.path.join(os.path.dirname(case_path), "common")
        has_case_specific_initial_stage = os.path.exists(
            os.path.join(case_path, "ProjectParameters_stage1.json")
        )
        stages = [1, 2] if has_case_specific_initial_stage else [2]
        KratosGeoMechanicsCPhiReductionProcess._run_stages_with_common_parameters(
            case_path, common_path, stages
        )
        reader = GiDOutputFileReader()
        actual_data = reader.read_output_from(os.path.join(case_path, "stage2.post.res"))
        node_ids = [57, 151, 235]
        return reader.nodal_values_at_time("DISPLACEMENT", 0.2, actual_data, node_ids)

    def _assert_case(self, case_name, expected_displacements_x):
        file_path = test_helper.get_file_path(
            os.path.join("C-Phi_reduction_process", case_name)
        )
        displacement_vectors = self._get_displacement_vectors(file_path)
        for index, expected_x in enumerate(expected_displacements_x):
            self.assertAlmostEqual(expected_x, displacement_vectors[index][0])

    def test_c_phi_reduction_process(self):
        test_cases = {
            "Mohr_Coulomb_model64": [-0.00198728, -0.00631146, -0.00450714],
            "Mohr_Coulomb_model": [-0.00222039, -0.00572449, -0.00338423],
            "Mohr_Coulomb_numerical_tangent": [-0.000371072, -0.000576557, 0.00110469],
            "Mohr_Coulomb_CL_APP": [-0.000371072, -0.000576557, 0.00110469]

        }

        for case_name, expected_displacements_x in test_cases.items():
            with self.subTest(case=case_name):
                self._assert_case(case_name, expected_displacements_x)

if __name__ == '__main__':
    KratosUnittest.main()
