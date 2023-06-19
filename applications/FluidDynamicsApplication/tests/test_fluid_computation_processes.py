import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as UnitTest
from KratosMultiphysics.testing.utilities import ReadModelPart
from KratosMultiphysics.process_factory import KratosProcessFactory

class FluidComputationProcessesTest(UnitTest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = Kratos.Model()
        cls.model_part = cls.model.CreateModelPart("test")
        cls.model_part.ProcessInfo[Kratos.DOMAIN_SIZE] = 2
        cls.model_part.AddNodalSolutionStepVariable(Kratos.VELOCITY)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.NORMAL)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.REACTION)
        ReadModelPart("AdjointVMSSensitivity2DTest/cylinder_test", cls.model_part)

        cls.density = 2.0
        cls.kinematic_viscosity = 3.0

        for element in cls.model_part.Elements:
            element.Properties[Kratos.DENSITY] = cls.density
            element.Properties[Kratos.DYNAMIC_VISCOSITY] = cls.kinematic_viscosity * cls.density
            break

        tmoc = Kratos.TetrahedralMeshOrientationCheck
        flags = (tmoc.COMPUTE_NODAL_NORMALS).AsFalse() | (tmoc.COMPUTE_CONDITION_NORMALS).AsFalse() | tmoc.ASSIGN_NEIGHBOUR_ELEMENTS_TO_CONDITIONS
        Kratos.TetrahedralMeshOrientationCheck(cls.model_part, False, flags).Execute()

        # now populate the velocity
        for node in cls.model_part.Nodes:
            node.SetSolutionStepValue(Kratos.VELOCITY, Kratos.Array3([node.Id, node.Id * 2, node.Id * 3]))
            node.SetSolutionStepValue(Kratos.REACTION, Kratos.Array3([node.Id, node.Id * 2, node.Id * 3]))

    def testComputeCFLProcess(self):
        settings = Kratos.Parameters(r'''
        [
            {
                "kratos_module" : "KratosMultiphysics.FluidDynamicsApplication",
                "python_module" : "compute_cfl_process",
                "Parameters" : {
                    "model_part_name" : "test"
                }
            }
        ]''')

        self.__RunTest(settings)

    def testComputeYPlusProcess(self):
        settings = Kratos.Parameters(r'''
        [
            {
                "kratos_module" : "KratosMultiphysics.FluidDynamicsApplication",
                "python_module" : "compute_y_plus_process",
                "Parameters" : {
                    "model_part_name" : "test"
                }
            }
        ]''')

        self.__RunTest(settings)

    def __RunTest(self, settings):
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
                process.Execute()                
                process.ExecuteFinalizeSolutionStep()

        for process in process_list:
            process.ExecuteFinalize()

if __name__ == '__main__':
    UnitTest.main()
