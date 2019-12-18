import KratosMultiphysics as Kratos
from KratosMultiphysics.process_factory import KratosProcessFactory

import KratosMultiphysics.KratosUnittest as UnitTest
import random


class TimeAveragingProcessTests(UnitTest.TestCase):
    def testTimeAveragingProcess(self):
        self.__CreateModel()

        settings = Kratos.Parameters(r'''
        [
            {
                "kratos_module" : "KratosMultiphysics",
                "python_module" : "process_factory",
                "process_name"  : "TimeAveragingProcess",
                "Parameters" : {
                    "model_part_name"                               : "test",
                    "variables_list"                                : ["VELOCITY", "DENSITY"],
                    "integration_start_point_control_variable_name" : "TIME",
                    "integration_start_point_control_value"         : 2.0
                }
            }
        ]''')

        factory = KratosProcessFactory(self.model)
        self.process_list = factory.ConstructListOfProcesses(settings)

        for process in self.process_list:
            process.Check()
        for process in self.process_list:
            process.ExecuteInitialize()

        velocity_vector = []
        density_vector = []
        for _ in self.model_part.Nodes:
            velocity_vector.append([0.0, 0.0, 0.0])
            density_vector.append(0.0)

        total_time = 0.0
        current_time = 0.0
        for _ in range(1, 10):
            current_time += 0.5
            self.model_part.CloneTimeStep(current_time)
            self.model_part.ProcessInfo[Kratos.STEP] += 1

            for process in self.process_list:
                process.ExecuteInitializeSolutionStep()

            for node_index, node in enumerate(self.model_part.Nodes):
                velocity = node.GetSolutionStepValue(Kratos.VELOCITY)
                temp_velocity = Kratos.Array3()
                temp_velocity[0] = velocity[0] * current_time
                temp_velocity[1] = velocity[1] * current_time
                temp_velocity[2] = velocity[2] * current_time
                node.SetSolutionStepValue(Kratos.VELOCITY, 0, temp_velocity)

                density = node.GetSolutionStepValue(Kratos.DENSITY)
                density = density * current_time
                node.SetSolutionStepValue(Kratos.DENSITY, 0, density)

                if (self.model_part.ProcessInfo[Kratos.TIME] >= 2.0):
                    velocity_vector[node_index][0] += temp_velocity[0] * 0.5
                    velocity_vector[node_index][1] += temp_velocity[1] * 0.5
                    velocity_vector[node_index][2] += temp_velocity[2] * 0.5
                    density_vector[node_index] += density * 0.5

            if (self.model_part.ProcessInfo[Kratos.TIME] >= 2.0):
                total_time += 0.5

            for process in self.process_list:
                process.ExecuteFinalizeSolutionStep()

        for process in self.process_list:
            process.ExecuteFinalize()

        for node_index, node in enumerate(self.model_part.Nodes):
            averaged_velocity = node.GetValue(Kratos.VELOCITY)
            averaged_density = node.GetValue(Kratos.DENSITY)

            self.assertAlmostEqual(averaged_velocity[0],
                                   velocity_vector[node_index][0] / total_time,
                                   12)
            self.assertAlmostEqual(averaged_velocity[1],
                                   velocity_vector[node_index][1] / total_time,
                                   12)
            self.assertAlmostEqual(averaged_velocity[2],
                                   velocity_vector[node_index][2] / total_time,
                                   12)
            self.assertAlmostEqual(averaged_density,
                                   density_vector[node_index] / total_time, 12)

    def __CreateModel(self):
        self.model = Kratos.Model()
        self.model_part = self.model.CreateModelPart("test")
        self.model_part.AddNodalSolutionStepVariable(Kratos.VELOCITY)
        self.model_part.AddNodalSolutionStepVariable(Kratos.DENSITY)

        self.model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
        self.model_part.CreateNewNode(2, 1.0, 0.0, 0.0)
        self.model_part.CreateNewNode(3, 2.0, 0.0, 0.0)
        self.model_part.CreateNewNode(4, 2.0, 1.0, 0.0)
        self.model_part.CreateNewNode(5, 1.5, 1.0, 0.0)
        self.model_part.CreateNewNode(6, 1.0, 1.0, 0.0)

        for node in self.model_part.Nodes:
            vector = Kratos.Vector(3)
            vector[0] = random.random()
            vector[1] = random.random()
            vector[2] = random.random()
            if (node.SolutionStepsDataHas(Kratos.VELOCITY)):
                vector = node.SetSolutionStepValue(Kratos.VELOCITY, 0, vector)

            scalar = random.random()
            if (node.SolutionStepsDataHas(Kratos.DENSITY)):
                node.SetSolutionStepValue(Kratos.DENSITY, 0, scalar)


if __name__ == '__main__':
    UnitTest.main()