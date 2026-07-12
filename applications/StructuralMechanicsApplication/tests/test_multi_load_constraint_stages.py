import pathlib

import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.StructuralMechanicsApplication

from KratosMultiphysics.modelers.import_mdpa_modeler import ImportMDPAModeler
from KratosMultiphysics.StructuralMechanicsApplication.multi_load_constraint_preparation import MultiLoadConstraintPreparation
from KratosMultiphysics.StructuralMechanicsApplication.multi_load_constraint_analysis import MultiLoadConstraintAnalysis

class TestMultiLoadConstraintStages(KratosUnittest.TestCase):

    _PREPARATION_STAGE_NAME = "prepare_lhs_and_rhss"
    _FIRST_ANALYSIS_STAGE_NAME = "solve_fixity_set_1"
    _SECOND_ANALYSIS_STAGE_NAME = "solve_fixity_set_2"

    def test_preparation_stage(self):
        with KratosUnittest.WorkFolderScope("multi_load_constraint_test", __file__):
            _, _, preparation_data = self._RunPreparation()

            lhs = preparation_data["lhs"]
            rhss = preparation_data["rhss"]
            fixities = preparation_data["fixities"]
            dofset = preparation_data["strategy_data"].GetDofSet()

            node_1_x_index = self._FindDofIndex(
                dofset,
                node_id=1,
                variable=KratosMultiphysics.DISPLACEMENT_X,
            )
            node_2_x_index = self._FindDofIndex(
                dofset,
                node_id=2,
                variable=KratosMultiphysics.DISPLACEMENT_X,
            )

            expected_stiffness = 210000.0 * 350.0 / 1000.0

            self.assertEqual(lhs.Size1(), 6)
            self.assertEqual(lhs.Size2(), 6)
            self.assertEqual(len(dofset), 6)

            self.assertAlmostEqual(lhs[node_1_x_index, node_1_x_index], expected_stiffness, places=12)
            self.assertAlmostEqual(lhs[node_1_x_index, node_2_x_index], -expected_stiffness, places=12)
            self.assertAlmostEqual(lhs[node_2_x_index, node_1_x_index], -expected_stiffness, places=12)
            self.assertAlmostEqual(lhs[node_2_x_index, node_2_x_index], expected_stiffness, places=12)

            self.assertIn("lc1", rhss)
            self.assertIn("lc2", rhss)
            self.assertAlmostEqual(rhss["lc1"][node_1_x_index], 1.0, places=12)
            self.assertAlmostEqual(rhss["lc1"][node_2_x_index], 0.0, places=12)
            self.assertAlmostEqual(rhss["lc2"][node_1_x_index], 0.0, places=12)
            self.assertAlmostEqual(rhss["lc2"][node_2_x_index], 1.0, places=12)

            self.assertEqual(len(fixities), 4)
            self._CheckFixityDefinition(
                fixities["fix1"],
                model_part_name="Structure.SPC_Group_Node1",
                constrained=[True, True, True],
            )
            self._CheckFixityDefinition(
                fixities["fix2"],
                model_part_name="Structure.SPC_Group_Node2",
                constrained=[False, True, True],
            )
            self._CheckFixityDefinition(
                fixities["fix3"],
                model_part_name="Structure.SPC_Group_Node1",
                constrained=[False, True, True],
            )
            self._CheckFixityDefinition(
                fixities["fix4"],
                model_part_name="Structure.SPC_Group_Node2",
                constrained=[True, True, True],
            )

            for dof in dofset:
                self.assertFalse(
                    dof.IsFixed(),
                    msg=(
                        "The preparation stage fixed "
                        f"{dof.GetVariable().Name()} on node {dof.Id()}."
                    ),
                )

    def test_two_consecutive_analysis_stages(self):
        with KratosUnittest.WorkFolderScope("multi_load_constraint_test", __file__):
            project_parameters, model, preparation_data = self._RunPreparation()

            first_solution = self._RunAnalysisStage(model, 
                                                    project_parameters, 
                                                    preparation_data, 
                                                    self._FIRST_ANALYSIS_STAGE_NAME)

            second_solution = self._RunAnalysisStage(model, 
                                                     project_parameters, 
                                                     preparation_data, 
                                                     self._SECOND_ANALYSIS_STAGE_NAME
            )

            dofset = preparation_data["strategy_data"].GetDofSet()

            node_1_x_index = self._FindDofIndex(dofset, 
                                                node_id=1, 
                                                variable=KratosMultiphysics.DISPLACEMENT_X,
            )
            node_2_x_index = self._FindDofIndex(dofset, 
                                                node_id=2, 
                                                variable=KratosMultiphysics.DISPLACEMENT_X,
            )

            young_modulus = 210000.0
            cross_area = 350.0
            length = 1000.0
            load = 10.0

            expected_displacement = (load * length / (young_modulus * cross_area))
            first_combination = first_solution[12]
            second_combination = second_solution[4]

            # First state:
            # node 1 fixed in x, node 2 free and loaded in x.
            self.assertAlmostEqual(
                first_combination[node_1_x_index],
                0.0,
                places=12,
            )
            self.assertAlmostEqual(
                first_combination[node_2_x_index],
                expected_displacement,
                places=12,
            )

            # Second state:
            # node 1 free and loaded in x, node 2 fixed in x.
            self.assertAlmostEqual(
                second_combination[node_1_x_index],
                expected_displacement,
                places=12,
            )
            self.assertAlmostEqual(
                second_combination[node_2_x_index],
                0.0,
                places=12,
            )

            # No stage-specific fixities should remain after both stages.
            for dof in dofset:
                self.assertFalse(
                    dof.IsFixed(),
                    msg=(
                        f"DOF {dof.GetVariable().Name()} on node "
                        f"{dof.Id()} remained fixed."
                    ),
                )

    def _RunPreparation(self):
        project_parameters = self._ReadParameters("ProjectParameters.json")
        model = KratosMultiphysics.Model()

        preparation_parameters = self._GetStageParameters(
            project_parameters,
            self._PREPARATION_STAGE_NAME,
        )
        preparation = MultiLoadConstraintPreparation(model, preparation_parameters)

        # The stage must add its nodal variables before nodes are imported.
        self._ImportModelPart(model)
        preparation.Run()

        return project_parameters, model, preparation.GetFinalData()

    def _CheckFixityDefinition(self, definition, model_part_name, constrained):
        parameters = definition["Parameters"]

        self.assertEqual(parameters["model_part_name"].GetString(), model_part_name)
        self.assertEqual(parameters["variable_name"].GetString(), "DISPLACEMENT")
        self.assertEqual(
            [value.GetBool() for value in parameters["constrained"].values()],
            constrained,
        )

    def _RunAnalysisStage(
        self,
        model,
        project_parameters,
        preparation_data,
        stage_name,
    ):
        analysis_parameters = self._GetStageParameters(
            project_parameters,
            stage_name,
        )

        analysis = MultiLoadConstraintAnalysis(
            model,
            analysis_parameters,
        )

        # Reproduce the data transfer normally performed by the orchestrator.
        output_data = {
            self._PREPARATION_STAGE_NAME: preparation_data
        }

        analysis.GetProjectOutputData(
            self._PREPARATION_STAGE_NAME,
            output_data,
        )
        analysis.Run()

        return analysis.GetFinalData()

    def _ImportModelPart(self, model):
        modeler_parameters = KratosMultiphysics.Parameters(r"""
        {
            "input_filename": "two_node_truss",
            "model_part_name": "Structure"
        }
        """)

        modeler = ImportMDPAModeler(model, modeler_parameters)
        modeler.SetupGeometryModel()
        modeler.PrepareGeometryModel()
        modeler.SetupModelPart()

    @staticmethod
    def _GetStageParameters(project_parameters, stage_name):
        return project_parameters[
            "stages"
        ][
            stage_name
        ][
            "stage_settings"
        ].Clone()

    @staticmethod
    def _ReadParameters(file_name):
        file_path = pathlib.Path(file_name)

        with file_path.open("r") as parameter_file:
            return KratosMultiphysics.Parameters(
                parameter_file.read()
            )

    def _FindDofIndex(self, dofset, node_id, variable):
        for index, dof in enumerate(dofset):
            if (
                dof.Id() == node_id
                and dof.GetVariable() == variable
            ):
                return index

        self.fail(
            f"Could not find {variable.Name()} on node {node_id}."
        )


if __name__ == "__main__":
    KratosUnittest.main()
