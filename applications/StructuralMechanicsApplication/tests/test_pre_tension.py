# --- External Imports ---
import numpy

# --- Kratos Imports ---
import KratosMultiphysics
import KratosMultiphysics.KratosUnittest
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import StructuralMechanicsAnalysis
from KratosMultiphysics.KratosUnittest import WorkFolderScope

# --- STD Imports ---
import pathlib
import json


class TestPretension(KratosMultiphysics.KratosUnittest.TestCase):

    def test_NeumannPreTensionedBeam(self) -> None:
        with WorkFolderScope("constraints", pathlib.Path(__file__).absolute()):
            # Load the configuration.
            with open("pre_tensioned_beam.json", "r") as file:
                configuration = json.loads(file.read())

            # Insert pre-tensioning.
            pre_tensioning_magnitude: float = 2e7
            configuration["processes"]["constraints_process_list"].append({
                "python_module" : "neumann_pre_tension_process",
                "kratos_module" : "KratosMultiphysics.StructuralMechanicsApplication",
                "Parameters" : {
                    "model_part_name" : "root.pre_tension_surface",
                    "magnitude" : pre_tensioning_magnitude,
                    "verbosity" : 3
                }})

            # Construct and run the analysis.
            parameters = KratosMultiphysics.Parameters(json.dumps(configuration))
            model = KratosMultiphysics.Model()
            analysis = StructuralMechanicsAnalysis(model, parameters)
            analysis.Run()

            # Check the reaction on the virtual nodes.
            model_part: KratosMultiphysics.ModelPart = model.GetModelPart("root").GetSubModelPart("pre_tension_surface")
            positive_side_virtual_node: KratosMultiphysics.Node = model_part.GetNode(8)
            negative_side_virtual_node: KratosMultiphysics.Node = model_part.GetNode(10)

            self.assertAlmostEqual(
                positive_side_virtual_node.GetSolutionStepValue(KratosMultiphysics.REACTION_X),
                -pre_tensioning_magnitude)
            self.assertAlmostEqual(
                negative_side_virtual_node.GetSolutionStepValue(KratosMultiphysics.REACTION_X),
                pre_tensioning_magnitude)


    def test_DirichletPreTensionedBeam(self) -> None:
        with WorkFolderScope("constraints", pathlib.Path(__file__).absolute()):
            # Load the configuration.
            with open("pre_tensioned_beam.json", "r") as file:
                configuration = json.loads(file.read())

            # Insert pre-tensioning.
            pre_tensioning_magnitude: float = 1e-1
            configuration["processes"]["constraints_process_list"].append({
                "python_module" : "dirichlet_pre_tension_process",
                "kratos_module" : "KratosMultiphysics.StructuralMechanicsApplication",
                "Parameters" : {
                    "model_part_name" : "root.pre_tension_surface",
                    "magnitude" : pre_tensioning_magnitude,
                    "verbosity" : 3
                }})

            # Construct and run the analysis.
            parameters = KratosMultiphysics.Parameters(json.dumps(configuration))
            model = KratosMultiphysics.Model()
            analysis = StructuralMechanicsAnalysis(model, parameters)
            analysis.Run()

            # Check the displacements on the duplicate and virtual nodes.
            model_part: KratosMultiphysics.ModelPart = model.GetModelPart("root")
            original_node: KratosMultiphysics.Node = model_part.GetNode(7)
            duplicate_node: KratosMultiphysics.Node = model_part.GetNode(9)

            normal: numpy.ndarray = numpy.array([2.0, 3.0, 0.0])
            normal /= numpy.linalg.norm(normal)

            relative_displacement: numpy.ndarray = numpy.array(
                duplicate_node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT)
                - original_node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT))
            relative_out_of_plane_displacement = numpy.dot(relative_displacement, normal) * normal
            relative_in_plane_displacement = relative_displacement - relative_out_of_plane_displacement

            #raise RuntimeError(relative_displacement)

            self.assertLessEqual(
                numpy.linalg.norm(relative_in_plane_displacement),
                1e-6)
            self.assertAlmostEqual(
                numpy.dot(relative_out_of_plane_displacement, normal),
                pre_tensioning_magnitude)


if __name__ == "__main__":
    KratosMultiphysics.KratosUnittest.main()
