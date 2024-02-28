import math
import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.testing.utilities import ReadModelPart
from KratosMultiphysics.OptimizationApplication.filtering.filtering_factory import Factory

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as kratos_unittest

class TestExplicitFilter(kratos_unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.model = Kratos.Model()
        cls.model_part = cls.model.CreateModelPart("test")
        cls.model_part.ProcessInfo[Kratos.DOMAIN_SIZE] = 3
        cls.model_part.AddNodalSolutionStepVariable(Kratos.NORMAL)
        ReadModelPart("shell", cls.model_part)

        Kratos.NormalCalculationUtils().CalculateNormals(cls.model_part, False, Kratos.NORMAL )

        for node in cls.model_part.Nodes:
            normal = node.GetSolutionStepValue(Kratos.NORMAL)
            normal_magnitude = (normal[0] ** 2 + normal[1] ** 2 + normal[2] ** 2) ** (0.5)
            node.SetValue(Kratos.VELOCITY, normal * (1 / normal_magnitude ** 2))

    def test_Hello(self):
        nodal_expression = Kratos.Expression.NodalExpression(self.model_part)
        Kratos.Expression.VariableExpressionIO.Read(nodal_expression, Kratos.VELOCITY, False)

        filtering_settings = Kratos.Parameters("""{
            "filter_type"               : "explicit_vertex_morphing",
            "filtering_model_part_name" : "test",
            "filter_function_type"      : "cosine",
            "radius"                    : 0.5,
            "max_nodes_in_filter_radius": 1000,
            "damping_function_type"     : "cosine",
            "damped_model_part_names"    : {
                //"DISPLACEMENT_edges" : [true, false, true]
                "test.DISPLACEMENT_points": [true, true, true]
            }
        }""")
        explicit_filter = Factory(self.model, "test", None, Kratos.Globals.DataLocation.NodeHistorical, filtering_settings)
        explicit_filter.Initialize()
        filtered_nodal_expression = explicit_filter.FilterField(nodal_expression)

        vtu_output = Kratos.VtuOutput(self.model_part)
        vtu_output.AddContainerExpression("raw", nodal_expression.Clone())
        vtu_output.AddContainerExpression("filtered", filtered_nodal_expression.Clone())
        vtu_output.PrintOutput("reza")

if __name__ == "__main__":
    kratos_unittest.main()