import KratosMultiphysics as Kratos
from KratosMultiphysics.process_factory import KratosProcessFactory

import KratosMultiphysics.StatisticsApplication as KratosStats
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.StatisticsApplication.test_utilities import HistoricalRetrievalMethod
from KratosMultiphysics.StatisticsApplication.test_utilities import NonHistoricalRetrievalMethod
from KratosMultiphysics.StatisticsApplication.test_utilities import InitializeContainerArrays
from KratosMultiphysics.StatisticsApplication.test_utilities import CheckValues
from KratosMultiphysics.StatisticsApplication.test_utilities import CreateModelPart
from KratosMultiphysics.StatisticsApplication.test_utilities import InitializeModelPartVariables
from KratosMultiphysics.StatisticsApplication.test_utilities import GetInitialVariableValue
from KratosMultiphysics.StatisticsApplication.test_utilities import InitializeProcesses
from KratosMultiphysics.StatisticsApplication.test_utilities import ExecuteProcessFinalizeSolutionStep


class TemporalSumMethodTests(KratosUnittest.TestCase):
    def setUp(self):
        self.model = Kratos.Model()
        self.model_part = self.model.CreateModelPart("test_model_part")
        self.model_part.SetBufferSize(1)

        self.__AddNodalSolutionStepVariables()
        CreateModelPart(self.model_part)
        InitializeModelPartVariables(self.model_part)

    def testSumHistoricalHistoricalValueMethod(self):
        norm_type = "none"
        settings = TemporalSumMethodTests.__GetDefaultSettings(
            norm_type, "nodal_historical_historical")
        self.__TestMethod(norm_type, settings, self.model_part.Nodes,
                          HistoricalRetrievalMethod, HistoricalRetrievalMethod)

    def testSumHistoricalNonHistoricalValueMethod(self):
        norm_type = "none"
        settings = TemporalSumMethodTests.__GetDefaultSettings(
            norm_type, "nodal_historical_non_historical")
        self.__TestMethod(norm_type, settings, self.model_part.Nodes,
                          HistoricalRetrievalMethod,
                          NonHistoricalRetrievalMethod)

    def testSumHistoricalHistoricalNormMethod(self):
        norm_type = "magnitude"
        settings = TemporalSumMethodTests.__GetDefaultSettings(
            norm_type, "nodal_historical_historical")
        self.__TestMethod(norm_type, settings, self.model_part.Nodes,
                          HistoricalRetrievalMethod, HistoricalRetrievalMethod)

    def testSumHistoricalNonHistoricalNormMethod(self):
        norm_type = "magnitude"
        settings = TemporalSumMethodTests.__GetDefaultSettings(
            norm_type, "nodal_historical_non_historical")
        self.__TestMethod(norm_type, settings, self.model_part.Nodes,
                          HistoricalRetrievalMethod,
                          NonHistoricalRetrievalMethod)

    def testSumNodalNonHistoricalValueMethod(self):
        norm_type = "none"
        settings = TemporalSumMethodTests.__GetDefaultSettings(
            norm_type, "nodal_non_historical")
        self.__TestMethod(norm_type, settings, self.model_part.Nodes,
                          NonHistoricalRetrievalMethod,
                          NonHistoricalRetrievalMethod)

    def testSumNodalNonHistoricalNormMethod(self):
        norm_type = "magnitude"
        settings = TemporalSumMethodTests.__GetDefaultSettings(
            norm_type, "nodal_non_historical")
        self.__TestMethod(norm_type, settings, self.model_part.Nodes,
                          NonHistoricalRetrievalMethod,
                          NonHistoricalRetrievalMethod)

    def testSumConditionNonHistoricalValueMethod(self):
        norm_type = "none"
        settings = TemporalSumMethodTests.__GetDefaultSettings(
            norm_type, "condition_non_historical")
        self.__TestMethod(norm_type, settings, self.model_part.Conditions,
                          NonHistoricalRetrievalMethod,
                          NonHistoricalRetrievalMethod)

    def testSumConditionNonHistoricalNormMethod(self):
        norm_type = "magnitude"
        settings = TemporalSumMethodTests.__GetDefaultSettings(
            norm_type, "condition_non_historical")
        self.__TestMethod(norm_type, settings, self.model_part.Conditions,
                          NonHistoricalRetrievalMethod,
                          NonHistoricalRetrievalMethod)

    def testSumElementNonHistoricalValueMethod(self):
        norm_type = "none"
        settings = TemporalSumMethodTests.__GetDefaultSettings(
            norm_type, "element_non_historical")
        self.__TestMethod(norm_type, settings, self.model_part.Elements,
                          NonHistoricalRetrievalMethod,
                          NonHistoricalRetrievalMethod)

    def testSumElementNonHistoricalNormMethod(self):
        norm_type = "magnitude"
        settings = TemporalSumMethodTests.__GetDefaultSettings(
            norm_type, "element_non_historical")
        self.__TestMethod(norm_type, settings, self.model_part.Elements,
                          NonHistoricalRetrievalMethod,
                          NonHistoricalRetrievalMethod)

    def __TestMethod(self, norm_type, settings, container, input_method,
                     output_method):
        factory = KratosProcessFactory(self.model)
        self.process_list = factory.ConstructListOfProcesses(settings)
        InitializeProcesses(self)

        scalar_list, vec_3d_list, vec_list, mat_list = InitializeContainerArrays(
            container)

        for step in range(0, 12, 2):
            self.model_part.CloneTimeStep(step)
            InitializeModelPartVariables(self.model_part)
            ExecuteProcessFinalizeSolutionStep(self)

            for index, item in enumerate(container):
                current_scalar = input_method(item, Kratos.PRESSURE)
                current_vector_3d = input_method(item, Kratos.VELOCITY)
                current_vector = input_method(item, Kratos.LOAD_MESHES)
                current_matrix = input_method(
                    item, Kratos.GREEN_LAGRANGE_STRAIN_TENSOR)

                if (step >= 4):
                    scalar_list[index].append(current_scalar)
                    vec_3d_list[index].append(current_vector_3d)
                    vec_list[index].append(current_vector)
                    mat_list[index].append(current_matrix)

                analytical_method_scalar = TemporalSumMethodTests.__AnalyticalMethod(
                    norm_type, Kratos.PRESSURE, scalar_list[index])
                analytical_method_vec_3d = TemporalSumMethodTests.__AnalyticalMethod(
                    norm_type, Kratos.VELOCITY, vec_3d_list[index])
                analytical_method_vec = TemporalSumMethodTests.__AnalyticalMethod(
                    norm_type, Kratos.LOAD_MESHES, vec_list[index])
                analytical_method_mat = TemporalSumMethodTests.__AnalyticalMethod(
                    norm_type, Kratos.GREEN_LAGRANGE_STRAIN_TENSOR,
                    mat_list[index])

                if (norm_type == "none"):
                    method_scalar = output_method(item,
                                                  KratosStats.PRESSURE_MEAN)
                    method_vec_3d = output_method(item,
                                                  KratosStats.VELOCITY_MEAN)
                    method_vec = output_method(item,
                                               Kratos.MATERIAL_PARAMETERS)
                    method_mat = output_method(item,
                                               Kratos.CAUCHY_STRESS_TENSOR)
                else:
                    method_scalar = output_method(item,
                                                  KratosStats.PRESSURE_NORM)
                    method_vec_3d = output_method(item,
                                                  KratosStats.VELOCITY_NORM)
                    method_vec = output_method(item, Kratos.DENSITY)
                    method_mat = output_method(item, Kratos.VISCOSITY)

                CheckValues(self, analytical_method_scalar, method_scalar, 8)
                CheckValues(self, analytical_method_vec_3d, method_vec_3d, 8)
                CheckValues(self, analytical_method_vec, method_vec, 16)
                CheckValues(self, analytical_method_mat, method_mat, 8)

    @staticmethod
    def __AnalyticalMethod(norm_type, variable, value_array):
        if (norm_type == "none"):
            result = GetInitialVariableValue(variable, "none")
            for item in value_array:
                result += item
        else:
            result = 0.0
            norm_method = KratosStats.MethodUtilities.GetNormMethod(
                variable, norm_type)
            for item in value_array:
                result += norm_method(item)

        return result * 2.0

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
            settings_str = settings_str.replace(
                "<OUTPUT_VARIABLES>",
                r'"VELOCITY_MEAN", "PRESSURE_MEAN", "MATERIAL_PARAMETERS", "CAUCHY_STRESS_TENSOR"'
            )
        else:
            settings_str = settings_str.replace(
                "<OUTPUT_VARIABLES>",
                r'"VELOCITY_NORM", "PRESSURE_NORM", "DENSITY", "VISCOSITY"')

        return Kratos.Parameters(settings_str)

    def __AddNodalSolutionStepVariables(self):
        # input variables
        self.model_part.AddNodalSolutionStepVariable(Kratos.PRESSURE)
        self.model_part.AddNodalSolutionStepVariable(Kratos.VELOCITY)
        self.model_part.AddNodalSolutionStepVariable(Kratos.LOAD_MESHES)
        self.model_part.AddNodalSolutionStepVariable(
            Kratos.GREEN_LAGRANGE_STRAIN_TENSOR)

        # output variables for output_1
        self.model_part.AddNodalSolutionStepVariable(KratosStats.PRESSURE_MEAN)
        self.model_part.AddNodalSolutionStepVariable(KratosStats.VELOCITY_MEAN)
        self.model_part.AddNodalSolutionStepVariable(
            Kratos.MATERIAL_PARAMETERS)
        self.model_part.AddNodalSolutionStepVariable(
            Kratos.CAUCHY_STRESS_TENSOR)

        self.model_part.AddNodalSolutionStepVariable(KratosStats.PRESSURE_NORM)
        self.model_part.AddNodalSolutionStepVariable(KratosStats.VELOCITY_NORM)
        self.model_part.AddNodalSolutionStepVariable(Kratos.DENSITY)
        self.model_part.AddNodalSolutionStepVariable(Kratos.VISCOSITY)


if __name__ == '__main__':
    KratosUnittest.main()