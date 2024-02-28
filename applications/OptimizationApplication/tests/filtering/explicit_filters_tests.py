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

    def test_Hello(self):
        filtering_settings = Kratos.Parameters("""{
            "filter_type"               : "explicit_vertex_morphing",
            "filtering_model_part_name" : "test",
            "filter_function_type"      : "linear",
            "radius"                    : 0.5,
            "max_nodes_in_filter_radius": 1000,
            "damping_function_type"     : "cosine",
            "damped_model_part_names"    : {
                //"test.DISPLACEMENT_edges" : [true, true, true]
                "test.DISPLACEMENT_points": [true, true, true]
            }
        }""")
        explicit_filter = Factory(self.model, "test", None, Kratos.Globals.DataLocation.NodeHistorical, filtering_settings)
        explicit_filter.Initialize()

        step_size = 5e-2
        for i in range(10):
            Kratos.NormalCalculationUtils().CalculateNormals(self.model_part, False, Kratos.NORMAL )
            nodal_expression = Kratos.Expression.NodalExpression(self.model_part)
            Kratos.Expression.VariableExpressionIO.Read(nodal_expression, Kratos.NORMAL, True)

            filtered_nodal_expression = explicit_filter.FilterField(nodal_expression)
            normalized_filtered = filtered_nodal_expression / Kratos.Expression.Utils.NormInf(filtered_nodal_expression)

            update = explicit_filter.FilterField(normalized_filtered * step_size)

            nodal_locs = Kratos.Expression.NodalExpression(self.model_part)
            Kratos.Expression.NodalPositionExpressionIO.Read(nodal_locs, Kratos.Configuration.Initial)
            final_pos = nodal_locs + update
            Kratos.Expression.NodalPositionExpressionIO.Write(final_pos, Kratos.Configuration.Initial)
            Kratos.Expression.NodalPositionExpressionIO.Write(final_pos, Kratos.Configuration.Current)

            explicit_filter.Update()

            vtu_output = Kratos.VtuOutput(self.model_part, binary_output=Kratos.VtuOutput.ASCII)
            vtu_output.AddContainerExpression("raw", nodal_expression.Clone())
            vtu_output.AddContainerExpression("filtered", filtered_nodal_expression.Clone())
            vtu_output.AddContainerExpression("normalized_filtered", normalized_filtered.Clone())
            vtu_output.AddContainerExpression("update", update.Clone())
            vtu_output.PrintOutput(f"output_{i}")

if __name__ == "__main__":
    kratos_unittest.main()
