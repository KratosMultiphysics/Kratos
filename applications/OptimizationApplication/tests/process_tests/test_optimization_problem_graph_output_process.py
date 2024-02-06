from pathlib import Path
import KratosMultiphysics as Kratos

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as kratos_unittest
import KratosMultiphysics.kratos_utilities as kratos_utils
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
try:
    from KratosMultiphysics.OptimizationApplication.processes.optimization_problem_graph_output_process import OptimizationProblemGraphOutputProcess
except:
    OptimizationProblemGraphOutputProcess = None

class TestOptimizationProblemGraphOutputProcess(kratos_unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.optimization_problem = OptimizationProblem(0)

        params = Kratos.Parameters("""
        {
            "output_file_name": "output_<step>",
            "max_iterations"  : 10,
            "output_interval": 5,
            "interval"       : [0, "End"],
            "x_axis_label"    : "Iteration [-]",
            "graphs"          : [
                {
                    "y_axis_label"   : "objs",
                    "y_min"          : "0.0",
                    "y_max"          : "<initial_max_value> * 20",
                    "is_log"         : false,
                    "legend_position": "upper right",
                    "components"     : ["test_component_1.value", "test_component_2.value"]
                },
                {
                    "y_axis_label"   : "constraints",
                    "y_min"          : "1.0",
                    "y_max"          : "<initial_max_value> * 20",
                    "is_log"         : true,
                    "legend_position": "upper right",
                    "components"     : ["test_component_3.value"]
                }
            ]
        }""")

        if OptimizationProblemGraphOutputProcess:
            cls.output_process = OptimizationProblemGraphOutputProcess(params, cls.optimization_problem)

        ComponentDataView("test_component_1", cls.optimization_problem).SetDataBuffer(1)
        ComponentDataView("test_component_2", cls.optimization_problem).SetDataBuffer(1)
        ComponentDataView("test_component_3", cls.optimization_problem).SetDataBuffer(1)

    @kratos_unittest.skipIf(OptimizationProblemGraphOutputProcess == None, "matplotlib is not found")
    def test_OptimizationProblemGraphOutputProcess(self):
        file_1 = Path("output_5.png")
        file_2 = Path("output_10.png")
        self.addCleanup(kratos_utils.DeleteFileIfExisting, str(file_1))
        self.addCleanup(kratos_utils.DeleteFileIfExisting, str(file_2))
        self.output_process.ExecuteInitialize()

        for step in range(10):
            self.optimization_problem.AdvanceStep()
            comp_data = ComponentDataView("test_component_1", self.optimization_problem).GetBufferedData()
            comp_data["value"] = step * 2 + 1
            comp_data = ComponentDataView("test_component_2", self.optimization_problem).GetBufferedData()
            comp_data["value"] = step * 3 + 1
            comp_data = ComponentDataView("test_component_3", self.optimization_problem).GetBufferedData()
            comp_data["value"] = step * 4 + 1
            self.output_process.ExecuteInitializeSolutionStep()
            self.output_process.ExecuteFinalizeSolutionStep()

            if self.output_process.IsOutputStep():
                self.output_process.PrintOutput()

        self.output_process.ExecuteFinalize()

        self.assertTrue(file_1.is_file())
        self.assertTrue(file_2.is_file())

if __name__ == "__main__":
    Kratos.Tester.SetVerbosity(Kratos.Tester.Verbosity.PROGRESS)  # TESTS_OUTPUTS
    kratos_unittest.main()