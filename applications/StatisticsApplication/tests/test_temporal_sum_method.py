import KratosMultiphysics as Kratos

import KratosMultiphysics as Kratos
from KratosMultiphysics.process_factory import KratosProcessFactory

import KratosMultiphysics.StatisticsApplication as KratosStats
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.StatisticsApplication.spatial_utilities import GetItemContainer
from KratosMultiphysics.StatisticsApplication.method_utilities import GetNormTypeContainer
from KratosMultiphysics.StatisticsApplication.method_utilities import GetMethod

from random import uniform


class TemporalMethodTests(KratosUnittest.TestCase):
    def setUp(self):
        self.model = Kratos.Model()
        self.model_part = self.model.CreateModelPart("test_model_part")
        self.model_part.SetBufferSize(10)

        TemporalMethodTests.__AddNodalSolutionStepVariables(self.model_part)
        TemporalMethodTests.__CreateModelPart(self.model_part)
        TemporalMethodTests.__InitializeVariables(self.model_part)

    def tearDown(self):
        # Code here will be placed AFTER every test in this TestCase.
        pass

    def testSumMethod(self):
        settings = Kratos.Parameters(r'''
        [
            {
                "kratos_module" : "KratosMultiphysics.StatisticsApplication",
                "python_module" : "temporal_statistics_process",
                "Parameters" : {
                    "model_part_name"                : "test_model_part",
                    "input_variable_settings" : [
                        {
                             "method_name"     : "sum",
                             "norm_type"       : "none",
                             "container"       : "nodal_historical_historical",
                             "method_settings" : {
                                 "input_variables"  : ["VELOCITY", "PRESSURE"],
                                 "output_variables" : ["VELOCITY_MEAN", "PRESSURE_MEAN"]
                             }
                        }
                    ],
                    "statistics_start_point_control_variable_name" : "TIME",
                    "statistics_start_point_control_value"         : 4.0
                }
            }
        ]''')

        factory = KratosProcessFactory(self.model)
        self.process_list = factory.ConstructListOfProcesses(settings)
        self.__Initializerocesses()

        velocity_list = []
        pressure_list = []
        for _ in self.model_part.Nodes:
            velocity_list.append([])
            pressure_list.append([])

        def analytical_method(value_array, variable):
            result = TemporalMethodTests.__GetInitialValue(variable, "none")
            for item in value_array:
                result += item
            return result

        for step in range(0, 12, 2):
            TemporalMethodTests.__InitializeVariables(self.model_part)
            self.__ExecuteFinalizeSolutionStep()

            for index, node  in enumerate(self.model_part.Nodes):
                current_velocity = node.GetSolutionStepValue(Kratos.VELOCITY)
                current_pressure = node.GetSolutionStepValue(Kratos.PRESSURE)
                if (step > 4):
                    velocity_list[index].append(current_velocity)
                    pressure_list[index].append(current_pressure)
                analytical_velocity = analytical_method(velocity_list[index], Kratos.VELOCITY)
                analytical_pressure = analytical_method(pressure_list[index], Kratos.PRESSURE)

                self.__CheckValues(analytical_pressure, node.GetSolutionStepValue(KratosStats.PRESSURE_MEAN), 1e-12)
                self.__CheckValues(analytical_velocity, node.GetSolutionStepValue(KratosStats.VELOCITY_MEAN), 1e-12)

            self.model_part.CloneTimeStep(step)

    def __Initializerocesses(self):
        for process in self.process_list:
            process.Check()
        for process in self.process_list:
            process.ExecuteInitialize()

    def __ExecuteFinalizeSolutionStep(self):
        for process in self.process_list:
            process.ExecuteFinalizeSolutionStep()

    def __CheckValues(self, value_a, value_b, tolerance):
        if (isinstance(value_a, tuple)):
            for i in range(len(value_a)):
                self.__CheckValues(value_a[i], value_b[i], tolerance)
        else:
            if (isinstance(value_a, Kratos.Matrix)):
                self.assertMatrixAlmostEqual(value_a, value_b, tolerance)
            elif (isinstance(value_a, Kratos.Vector)):
                self.assertVectorAlmostEqual(value_a, value_b, tolerance)
            elif (isinstance(value_a, list)):
                for i in range(len(value_a)):
                    self.assertAlmostEqual(value_a[i], value_b[i], tolerance)
            else:
                self.assertAlmostEqual(value_a, value_b, tolerance)

    @staticmethod
    def __AddNodalSolutionStepVariables(model_part):
        model_part.AddNodalSolutionStepVariable(Kratos.PRESSURE)
        model_part.AddNodalSolutionStepVariable(Kratos.VELOCITY)
        model_part.AddNodalSolutionStepVariable(KratosStats.PRESSURE_MEAN)
        model_part.AddNodalSolutionStepVariable(KratosStats.VELOCITY_MEAN)
        model_part.AddNodalSolutionStepVariable(KratosStats.PRESSURE_VARIANCE)
        model_part.AddNodalSolutionStepVariable(KratosStats.VELOCITY_VARIANCE)
        model_part.AddNodalSolutionStepVariable(KratosStats.PRESSURE_NORM)
        model_part.AddNodalSolutionStepVariable(KratosStats.VELOCITY_NORM)
        model_part.AddNodalSolutionStepVariable(Kratos.LOAD_MESHES)
        model_part.AddNodalSolutionStepVariable(Kratos.GREEN_LAGRANGE_STRAIN_TENSOR)


    @staticmethod
    def __CreateModelPart(model_part):
        model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
        model_part.CreateNewNode(2, 1.0, 0.0, 0.0)
        model_part.CreateNewNode(3, 2.0, 0.0, 0.0)
        model_part.CreateNewNode(4, 2.0, 1.0, 0.0)
        model_part.CreateNewNode(5, 1.5, 1.0, 0.0)
        model_part.CreateNewNode(6, 0.5, 1.0, 0.0)
        model_part.CreateNewNode(7, 0.0, 1.0, 0.0)

        prop = model_part.GetProperties()[0]

        model_part.CreateNewElement("Element2D3N", 1, [1, 6, 7], prop)
        model_part.CreateNewElement("Element2D3N", 2, [1, 2, 6], prop)
        model_part.CreateNewElement("Element2D3N", 3, [6, 2, 5], prop)
        model_part.CreateNewElement("Element2D3N", 4, [2, 3, 5], prop)
        model_part.CreateNewElement("Element2D3N", 5, [5, 3, 4], prop)

        model_part.CreateNewCondition("Condition2D2N", 1, [1, 2], prop)
        model_part.CreateNewCondition("Condition2D2N", 2, [2, 3], prop)
        model_part.CreateNewCondition("Condition2D2N", 3, [3, 4], prop)
        model_part.CreateNewCondition("Condition2D2N", 4, [4, 5], prop)
        model_part.CreateNewCondition("Condition2D2N", 5, [5, 6], prop)
        model_part.CreateNewCondition("Condition2D2N", 6, [6, 7], prop)
        model_part.CreateNewCondition("Condition2D2N", 7, [7, 1], prop)

    @staticmethod
    def __InitializeVariables(model_part):

        for node in model_part.Nodes:
            s, a, v, m = TemporalMethodTests.__GetNewValue()
            node.SetValue(Kratos.PRESSURE, s)
            node.SetValue(Kratos.VELOCITY, a)
            node.SetValue(Kratos.LOAD_MESHES, v)
            node.SetValue(Kratos.GREEN_LAGRANGE_STRAIN_TENSOR, m)

            s, a, v, m = TemporalMethodTests.__GetNewValue()
            node.SetSolutionStepValue(Kratos.PRESSURE, 0, s)
            node.SetSolutionStepValue(Kratos.VELOCITY, 0, a)
            node.SetSolutionStepValue(Kratos.LOAD_MESHES, 0, v)
            node.SetSolutionStepValue(Kratos.GREEN_LAGRANGE_STRAIN_TENSOR, 0,
                                      m)

        for condition in model_part.Conditions:
            s, a, v, m = TemporalMethodTests.__GetNewValue()
            condition.SetValue(Kratos.PRESSURE, s)
            condition.SetValue(Kratos.VELOCITY, a)
            condition.SetValue(Kratos.LOAD_MESHES, v)
            condition.SetValue(Kratos.GREEN_LAGRANGE_STRAIN_TENSOR, m)

        for element in model_part.Elements:
            s, a, v, m = TemporalMethodTests.__GetNewValue()
            element.SetValue(Kratos.PRESSURE, s)
            element.SetValue(Kratos.VELOCITY, a)
            element.SetValue(Kratos.LOAD_MESHES, v)
            element.SetValue(Kratos.GREEN_LAGRANGE_STRAIN_TENSOR, m)

    @staticmethod
    def __GetValue(item, container_type, variable):
        if (container_type.endswith("non_historical")):
            return item.GetValue(variable)
        else:
            return item.GetSolutionStepValue(variable)

    @staticmethod
    def __GetNormValue(variable, value, norm_type):
        if (norm_type == "none"):
            return value

        norm_method = KratosStats.MethodUtilities.GetNormMethod(
            variable, norm_type)
        return norm_method(value)

    @staticmethod
    def __GetContainer(model_part, container_type):
        if (container_type.startswith("nodal")):
            return model_part.Nodes
        elif (container_type.startswith("element")):
            return model_part.Elements
        elif (container_type.startswith("condition")):
            return model_part.Conditions

    @staticmethod
    def __GetNewValue():
        m = Kratos.Matrix(5, 5)
        v = Kratos.Vector(5)
        a = Kratos.Vector(3)

        for i in range(5):
            v[i] = uniform(-10.0, 10.0)
            a[i % 3] = uniform(-10.0, 10.0)
            for j in range(5):
                m[i, j] = uniform(-10.0, 10.0)

        return uniform(-10.0, 10.0), a, v, m

    @staticmethod
    def __GetInitialValue(variable, norm_type):
        if (norm_type == "none"):
            if (variable == Kratos.PRESSURE):
                return 0.0
            elif (variable == Kratos.VELOCITY):
                return Kratos.Vector(3, 0.0)
            elif (variable == Kratos.LOAD_MESHES):
                return Kratos.Vector(5, 0.0)
            elif (variable == Kratos.GREEN_LAGRANGE_STRAIN_TENSOR):
                return Kratos.Matrix(5, 5, 0.0)
        else:
            return 0.0


if __name__ == '__main__':
    KratosUnittest.main()