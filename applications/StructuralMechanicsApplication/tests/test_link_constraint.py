# --- Kratos Imports ---
import KratosMultiphysics
import KratosMultiphysics.KratosUnittest
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import StructuralMechanicsAnalysis
from KratosMultiphysics.KratosUnittest import WorkFolderScope

# --- STD Imports ---
import pathlib
import math


class TestLinkConstraint(KratosMultiphysics.KratosUnittest.TestCase):

    def test_LinkConstraint2DWithoutMovedMesh(self) -> None:
        self.__Run2D(False)

    def test_LinkConstraint2DWithMovedMesh(self) -> None:
        self.__Run2D(True)

    def test_LinkConstraint3DWithoutMovedMesh(self) -> None:
        self.__Run3D(False)

    def test_LinkConstraint3DWithMovedMesh(self) -> None:
        self.__Run3D(True)

    def test_MixedConstraints(self) -> None:
        with WorkFolderScope("constraints", pathlib.Path(__file__).absolute()):
            with open("mixed_constraints.json") as project_parameters_file:
                parameters = KratosMultiphysics.Parameters(project_parameters_file.read())


            model = KratosMultiphysics.Model()
            analysis = StructuralMechanicsAnalysis(model, parameters)
            self.__ReadMDPA(model, "mixed_constraints")

            root_model_part = model.GetModelPart("root")
            link_constraint_node_ids = (2, 5)
            first: KratosMultiphysics.Node = root_model_part.GetNode(link_constraint_node_ids[0])
            second: KratosMultiphysics.Node = root_model_part.GetNode(link_constraint_node_ids[1])
            initial_distance: float = self.__GetNodeDistance(first, second, initial_configuration = True)

            analysis.Run()

            # Check whether the link constraint is satisfied.
            distance: float = self.__GetNodeDistance(first, second)
            self.assertAlmostEqual(distance, initial_distance, places = 3)

            # Check whether the master-slave constraints are satisfied.
            master_slave_constraint_node_ids = (3, 6)
            first = root_model_part.GetNode(master_slave_constraint_node_ids[0])
            second = root_model_part.GetNode(master_slave_constraint_node_ids[1])
            self.assertAlmostEqual(
                first.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X) / second.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X),
                2.0)
            self.assertAlmostEqual(
                first.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y) / second.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y),
                3.0)


    @staticmethod
    def __GetNodeDistance(first: KratosMultiphysics.Node,
                          second: KratosMultiphysics.Node,
                          initial_configuration: bool = False) -> float:
        initial_positions: "list[list[float]]" = [[first.X0, first.Y0, first.Z0],
                                                  [second.X0, second.Y0, second.Z0]]

        if not initial_configuration:
            displacements: "list[list[float]]" = [[first.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X),
                                                   first.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y),
                                                   first.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Z)],
                                                  [second.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X),
                                                   second.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y),
                                                   second.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Z)]]
        else:
            displacements: "list[list[float]]" = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]

        distance: float = 0.0
        for i_component in range(3):
            diff: float = initial_positions[0][i_component] - initial_positions[1][i_component]
            if not initial_configuration:
                diff += displacements[0][i_component] - displacements[1][i_component]
            distance += diff * diff
        return math.sqrt(distance)


    def __Run2D(self, move_mesh_flag: bool) -> None:
        with WorkFolderScope("constraints", pathlib.Path(__file__).absolute()):
            # Load config.
            with open("link_constraint_2d.json") as project_parameters_file:
                parameters = KratosMultiphysics.Parameters(project_parameters_file.read())

            parameters["solver_settings"]["move_mesh_flag"].SetBool(move_mesh_flag)
            parameters["processes"]["constraints_process_list"][1]["Parameters"]["move_mesh_flag"].SetBool(move_mesh_flag)

            # The test is set up such that 2 pseudo time steps are taken.
            #
            # At the end of the first one, the constraint is supposed to
            # be active so distance between the constrained nodes is supposed
            # to remain constant.
            #
            # At the end of the second step, the constraint is supposed to be
            # inactive, allowing the distance between the constrained nodes to
            # change.

            # Run the analysis until the first step and check whether constraints are satisfied.
            parameters["problem_data"]["end_time"].SetDouble(0.5)
            model = KratosMultiphysics.Model()
            analysis = StructuralMechanicsAnalysis(model, parameters)
            self.__ReadMDPA(model, "link_constraint_2d")

            root_model_part = model.GetModelPart("root")
            constrained_node_ids: "tuple[int,int]" = (2, 3)
            first: KratosMultiphysics.Node = root_model_part.GetNode(constrained_node_ids[0])
            second: KratosMultiphysics.Node = root_model_part.GetNode(constrained_node_ids[1])
            initial_distance: float = self.__GetNodeDistance(first, second, initial_configuration = True)

            analysis.Run()
            distance: float = self.__GetNodeDistance(first, second)
            self.assertAlmostEqual(distance, initial_distance, places = 3)

            # Run the analysis until the first step and check whether constraints are not satisfied anymore.
            parameters["problem_data"]["end_time"].SetDouble(1.0)
            model = KratosMultiphysics.Model()
            analysis = StructuralMechanicsAnalysis(model, parameters)
            self.__ReadMDPA(model, "link_constraint_2d")

            root_model_part = model.GetModelPart("root")
            constrained_node_ids: "tuple[int,int]" = (2, 3)
            first: KratosMultiphysics.Node = root_model_part.GetNode(constrained_node_ids[0])
            second: KratosMultiphysics.Node = root_model_part.GetNode(constrained_node_ids[1])
            initial_distance: float = self.__GetNodeDistance(first, second, initial_configuration = True)

            analysis.Run()
            distance: float = self.__GetNodeDistance(first, second)
            self.assertNotAlmostEqual(distance, initial_distance, places = 0)


    def __Run3D(self, move_mesh_flag: bool) -> None:
        with WorkFolderScope("constraints", pathlib.Path(__file__).absolute()):
            # Load config.
            with open("link_constraint_3d.json") as project_parameters_file:
                parameters = KratosMultiphysics.Parameters(project_parameters_file.read())

            parameters["solver_settings"]["move_mesh_flag"].SetBool(move_mesh_flag)
            parameters["processes"]["constraints_process_list"][1]["Parameters"]["move_mesh_flag"].SetBool(move_mesh_flag)

            model = KratosMultiphysics.Model()
            analysis = StructuralMechanicsAnalysis(model, parameters)
            self.__ReadMDPA(model, "link_constraint_3d")

            root_model_part = model.GetModelPart("root")
            constrained_node_ids: "tuple[int,int]" = (2, 3)
            first: KratosMultiphysics.Node = root_model_part.GetNode(constrained_node_ids[0])
            second: KratosMultiphysics.Node = root_model_part.GetNode(constrained_node_ids[1])
            initial_distance: float = self.__GetNodeDistance(first, second, initial_configuration = True)

            analysis.Run()
            distance: float = self.__GetNodeDistance(first, second)
            self.assertAlmostEqual(distance, initial_distance, places = 3)


    def __Run(self,
              analysis_stage: StructuralMechanicsAnalysis) -> None:
        analysis_stage.Run()


    @staticmethod
    def __ReadMDPA(model: KratosMultiphysics.Model,
                   mdpa_name: str) -> KratosMultiphysics.ModelPart:
        mdpa_path = pathlib.Path(mdpa_name)
        root_model_part = model.CreateModelPart("root")
        KratosMultiphysics.ModelPartIO(mdpa_path, KratosMultiphysics.ModelPartIO.READ).ReadModelPart(root_model_part)
        return root_model_part


if __name__ == "__main__":
    KratosMultiphysics.KratosUnittest.main()
