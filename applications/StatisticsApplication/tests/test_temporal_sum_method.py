import KratosMultiphysics as Kratos

import KratosMultiphysics as Kratos
from KratosMultiphysics.process_factory import KratosProcessFactory

import KratosMultiphysics.StatisticsApplication as KratosStats
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.StatisticsApplication.spatial_utilities import GetItemContainer
from KratosMultiphysics.StatisticsApplication.method_utilities import GetNormTypeContainer
from KratosMultiphysics.StatisticsApplication.method_utilities import GetMethod
from KratosMultiphysics.StatisticsApplication.test_utilities import HistoricalRetrievalMethod
from KratosMultiphysics.StatisticsApplication.test_utilities import NonHistoricalRetrievalMethod
from random import uniform


class TemporalSumMethodTests(KratosUnittest.TestCase):
    def setUp(self):
        self.model = Kratos.Model()
        self.model_part = self.model.CreateModelPart("test_model_part")
        self.model_part.SetBufferSize(10)

        TemporalSumMethodTests.__AddNodalSolutionStepVariables(self.model_part)
        TemporalSumMethodTests.__CreateModelPart(self.model_part)
        TemporalSumMethodTests.__InitializeVariables(self.model_part)

    def testSumHistoricalHistoricalValueMethod(self):
        norm_type = "none"
        settings = TemporalSumMethodTests.__GetDefaultSettings(norm_type, "nodal_historical_historical")
        self.__TestMethod(norm_type, settings, self.model_part.Nodes, HistoricalRetrievalMethod, HistoricalRetrievalMethod)

    def testSumHistoricalNonHistoricalValueMethod(self):
        norm_type = "none"
        settings = TemporalSumMethodTests.__GetDefaultSettings(norm_type, "nodal_historical_non_historical")
        self.__TestMethod(norm_type, settings, self.model_part.Nodes, HistoricalRetrievalMethod, NonHistoricalRetrievalMethod)

    def testSumHistoricalHistoricalNormMethod(self):
        norm_type = "magnitude"
        settings = TemporalSumMethodTests.__GetDefaultSettings(norm_type, "nodal_historical_historical")
        self.__TestMethod(norm_type, settings, self.model_part.Nodes, HistoricalRetrievalMethod, HistoricalRetrievalMethod)

    def testSumHistoricalNonHistoricalNormMethod(self):
        norm_type = "magnitude"
        settings = TemporalSumMethodTests.__GetDefaultSettings(norm_type, "nodal_historical_non_historical")
        self.__TestMethod(norm_type, settings, self.model_part.Nodes, HistoricalRetrievalMethod, NonHistoricalRetrievalMethod)

    def testSumNodalNonHistoricalValueMethod(self):
        norm_type = "none"
        settings = TemporalSumMethodTests.__GetDefaultSettings(norm_type, "nodal_non_historical")
        self.__TestMethod(norm_type, settings, self.model_part.Nodes, NonHistoricalRetrievalMethod, NonHistoricalRetrievalMethod)

    def testSumNodalNonHistoricalNormMethod(self):
        norm_type = "magnitude"
        settings = TemporalSumMethodTests.__GetDefaultSettings(norm_type, "nodal_non_historical")
        self.__TestMethod(norm_type, settings, self.model_part.Nodes, NonHistoricalRetrievalMethod, NonHistoricalRetrievalMethod)

    def testSumConditionNonHistoricalValueMethod(self):
        norm_type = "none"
        settings = TemporalSumMethodTests.__GetDefaultSettings(norm_type, "condition_non_historical")
        self.__TestMethod(norm_type, settings, self.model_part.Conditions, NonHistoricalRetrievalMethod, NonHistoricalRetrievalMethod)

    def testSumConditionNonHistoricalNormMethod(self):
        norm_type = "magnitude"
        settings = TemporalSumMethodTests.__GetDefaultSettings(norm_type, "condition_non_historical")
        self.__TestMethod(norm_type, settings, self.model_part.Conditions, NonHistoricalRetrievalMethod, NonHistoricalRetrievalMethod)

    def testSumElementNonHistoricalValueMethod(self):
        norm_type = "none"
        settings = TemporalSumMethodTests.__GetDefaultSettings(norm_type, "element_non_historical")
        self.__TestMethod(norm_type, settings, self.model_part.Elements, NonHistoricalRetrievalMethod, NonHistoricalRetrievalMethod)

    def testSumElementNonHistoricalNormMethod(self):
        norm_type = "magnitude"
        settings = TemporalSumMethodTests.__GetDefaultSettings(norm_type, "element_non_historical")
        self.__TestMethod(norm_type, settings, self.model_part.Elements, NonHistoricalRetrievalMethod, NonHistoricalRetrievalMethod)


    def __TestMethod(self, norm_type, settings, container, input_method, output_method):
        factory = KratosProcessFactory(self.model)
        self.process_list = factory.ConstructListOfProcesses(settings)
        self.__Initializerocesses()

        scalar_list, vec_3d_list, vec_list, mat_list = TemporalSumMethodTests.__InitializeArrays(container)

        for step in range(0, 12, 2):
            TemporalSumMethodTests.__InitializeVariables(self.model_part)
            self.__ExecuteFinalizeSolutionStep()

            for index, item in enumerate(container):
                current_scalar = input_method(item, Kratos.PRESSURE)
                current_vector_3d = input_method(item, Kratos.VELOCITY)
                current_vector = input_method(item, Kratos.LOAD_MESHES)
                current_matrix = input_method(item, Kratos.GREEN_LAGRANGE_STRAIN_TENSOR)

                if (step > 4):
                    scalar_list[index].append(current_scalar)
                    vec_3d_list[index].append(current_vector_3d)
                    vec_list[index].append(current_vector)
                    mat_list[index].append(current_matrix)

                analytical_method_scalar = TemporalSumMethodTests.__AnalyticalMethod(norm_type, Kratos.PRESSURE, scalar_list[index])
                analytical_method_vec_3d = TemporalSumMethodTests.__AnalyticalMethod(norm_type, Kratos.VELOCITY, vec_3d_list[index])
                analytical_method_vec = TemporalSumMethodTests.__AnalyticalMethod(norm_type, Kratos.LOAD_MESHES, vec_list[index])
                analytical_method_mat = TemporalSumMethodTests.__AnalyticalMethod(norm_type, Kratos.GREEN_LAGRANGE_STRAIN_TENSOR, mat_list[index])

                if (norm_type == "none"):
                    method_scalar = output_method(item, KratosStats.PRESSURE_MEAN)
                    method_vec_3d = output_method(item, KratosStats.VELOCITY_MEAN)
                    method_vec = output_method(item, Kratos.MATERIAL_PARAMETERS)
                    method_mat = output_method(item, Kratos.CAUCHY_STRESS_TENSOR)
                else:
                    method_scalar = output_method(item, KratosStats.PRESSURE_NORM)
                    method_vec_3d = output_method(item, KratosStats.VELOCITY_NORM)
                    method_vec = output_method(item, Kratos.DENSITY)
                    method_mat = output_method(item, Kratos.VISCOSITY)

                self.__CheckValues(analytical_method_scalar, method_scalar, 8)
                self.__CheckValues(analytical_method_vec_3d, method_vec_3d, 8)
                self.__CheckValues(analytical_method_vec, method_vec, 16)
                self.__CheckValues(analytical_method_mat, method_mat, 8)

            self.model_part.CloneTimeStep(step)

    @staticmethod
    def __InitializeArrays(container):
        scalar_list = []
        vector_3d_list = []
        vector_list = []
        matrix_list = []
        for _ in container:
            scalar_list.append([])
            vector_3d_list.append([])
            vector_list.append([])
            matrix_list.append([])

        return scalar_list, vector_3d_list, vector_list, matrix_list

    @staticmethod
    def __AnalyticalMethod(norm_type, variable, value_array):
        if (norm_type == "none"):
            result = TemporalSumMethodTests.__GetInitialValue(variable, "none")
            for item in value_array:
                result += item
        else:
            result = 0.0
            norm_method = KratosStats.MethodUtilities.GetNormMethod(variable, norm_type)
            for item in value_array:
                result += norm_method(item)

        return result

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
    def __GetDefaultSettings(norm_type, container_name):
        settings_str = r'''
        [
            {
                "kratos_module" : "KratosMultiphysics.StatisticsApplication",
                "python_module" : "temporal_statistics_process",
                "Parameters" : {
                    "model_part_name"                : "test_model_part",
                    "input_variable_settings" : [
                        {
                             "method_name"     : "sum",
                             "norm_type"       : "<TEST_NORM_TYPE>",
                             "container"       : "<TEST_CONTAINER>",
                             "echo_level"      : 0,
                             "method_settings" : {
                                 "input_variables"  : ["VELOCITY", "PRESSURE", "LOAD_MESHES", "GREEN_LAGRANGE_STRAIN_TENSOR"],
                                 "output_variables" : [<OUTPUT_VARIABLES>]
                             }
                        }
                    ],
                    "statistics_start_point_control_variable_name" : "TIME",
                    "statistics_start_point_control_value"         : 4.0
                }
            }
        ]'''
        settings_str = settings_str.replace("<TEST_NORM_TYPE>", norm_type)
        settings_str = settings_str.replace("<TEST_CONTAINER>", container_name)
        if (norm_type == "none"):
            settings_str = settings_str.replace("<OUTPUT_VARIABLES>", r'"VELOCITY_MEAN", "PRESSURE_MEAN", "MATERIAL_PARAMETERS", "CAUCHY_STRESS_TENSOR"')
        else:
            settings_str = settings_str.replace("<OUTPUT_VARIABLES>", r'"VELOCITY_NORM", "PRESSURE_NORM", "DENSITY", "VISCOSITY"')

        return Kratos.Parameters(settings_str)

    @staticmethod
    def __AddNodalSolutionStepVariables(model_part):
        # input variables
        model_part.AddNodalSolutionStepVariable(Kratos.PRESSURE)
        model_part.AddNodalSolutionStepVariable(Kratos.VELOCITY)
        model_part.AddNodalSolutionStepVariable(Kratos.LOAD_MESHES)
        model_part.AddNodalSolutionStepVariable(
            Kratos.GREEN_LAGRANGE_STRAIN_TENSOR)

        # output variables for output_1
        model_part.AddNodalSolutionStepVariable(KratosStats.PRESSURE_MEAN)
        model_part.AddNodalSolutionStepVariable(KratosStats.VELOCITY_MEAN)
        model_part.AddNodalSolutionStepVariable(Kratos.MATERIAL_PARAMETERS)
        model_part.AddNodalSolutionStepVariable(Kratos.CAUCHY_STRESS_TENSOR)

        model_part.AddNodalSolutionStepVariable(KratosStats.PRESSURE_NORM)
        model_part.AddNodalSolutionStepVariable(KratosStats.VELOCITY_NORM)
        model_part.AddNodalSolutionStepVariable(Kratos.DENSITY)
        model_part.AddNodalSolutionStepVariable(Kratos.VISCOSITY)

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
            s, a, v, m = TemporalSumMethodTests.__GetNewValue()
            node.SetValue(Kratos.PRESSURE, s)
            node.SetValue(Kratos.VELOCITY, a)
            node.SetValue(Kratos.LOAD_MESHES, v)
            node.SetValue(Kratos.GREEN_LAGRANGE_STRAIN_TENSOR, m)

            s, a, v, m = TemporalSumMethodTests.__GetNewValue()
            node.SetSolutionStepValue(Kratos.PRESSURE, 0, s)
            node.SetSolutionStepValue(Kratos.VELOCITY, 0, a)
            node.SetSolutionStepValue(Kratos.LOAD_MESHES, 0, v)
            node.SetSolutionStepValue(Kratos.GREEN_LAGRANGE_STRAIN_TENSOR, 0,
                                      m)

        for condition in model_part.Conditions:
            s, a, v, m = TemporalSumMethodTests.__GetNewValue()
            condition.SetValue(Kratos.PRESSURE, s)
            condition.SetValue(Kratos.VELOCITY, a)
            condition.SetValue(Kratos.LOAD_MESHES, v)
            condition.SetValue(Kratos.GREEN_LAGRANGE_STRAIN_TENSOR, m)

        for element in model_part.Elements:
            s, a, v, m = TemporalSumMethodTests.__GetNewValue()
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