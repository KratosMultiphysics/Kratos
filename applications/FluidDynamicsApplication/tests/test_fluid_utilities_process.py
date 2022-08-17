import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as UnitTest
from KratosMultiphysics.testing.utilities import ReadModelPart
from KratosMultiphysics.process_factory import KratosProcessFactory

class FluidUtilitiesProcessTest(UnitTest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = Kratos.Model()
        cls.model_part = cls.model.CreateModelPart("test")
        cls.model_part.ProcessInfo[Kratos.DOMAIN_SIZE] = 2
        cls.model_part.AddNodalSolutionStepVariable(Kratos.VELOCITY)
        ReadModelPart("AdjointVMSSensitivity2DTest/cylinder_test", cls.model_part)

        # now populate the velocity
        for node in cls.model_part.Nodes:
            node.SetSolutionStepValue(Kratos.VELOCITY, Kratos.Array3([node.Id, node.Id * 2, node.Id * 3]))

    def testCFLUtility(self):
        settings = Kratos.Parameters(r'''
        [
            {
                "kratos_module" : "KratosMultiphysics.FluidDynamicsApplication",
                "python_module" : "fluid_utilities_process",
                "Parameters" : {
                    "model_part_name" : "test",
                    "interval"        : [0.0, "End"],
                    "utility_settings": {
                        "utility_type" : "cfl",
                        "echo_level"   : 0
                    }
                }
            }
        ]''')

        factory = KratosProcessFactory(self.model)
        process_list = factory.ConstructListOfProcesses(settings)

        for process in process_list:
            process.Check()

        time_steps = [1.0, 1.5, 2.0, 2.5]
        for step, time_step in enumerate(time_steps, 1):
            self.model_part.ProcessInfo[Kratos.TIME] = time_step
            self.model_part.ProcessInfo[Kratos.STEP] = step

            for process in process_list:
                process.ExecuteInitializeSolutionStep()
                process.ExecuteFinalizeSolutionStep()

        for process in process_list:
            process.ExecuteFinalize()

if __name__ == '__main__':
    UnitTest.main()
