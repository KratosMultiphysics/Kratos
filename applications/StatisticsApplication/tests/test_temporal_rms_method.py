import KratosMultiphysics as Kratos
from KratosMultiphysics.process_factory import KratosProcessFactory

import KratosMultiphysics.StatisticsApplication as KratosStats
from KratosMultiphysics.StatisticsApplication.test_utilities import InitializeContainerArrays
from KratosMultiphysics.StatisticsApplication.test_utilities import CheckValues
from KratosMultiphysics.StatisticsApplication.test_utilities import GetInitialVariableValue
from KratosMultiphysics.StatisticsApplication.test_utilities import InitializeModelPartVariables
from KratosMultiphysics.StatisticsApplication.test_utilities import InitializeProcesses
from KratosMultiphysics.StatisticsApplication.test_utilities import ExecuteProcessFinalizeSolutionStep

import temporal_statistics_test_case


class TemporalRootMeanSquareMethodHelperClass(
        temporal_statistics_test_case.TemporalStatisticsTestCase):
    def RunTemporalStatisticsTest(self, norm_type, container_name):

        settings = TemporalRootMeanSquareMethodHelperClass.__GetDefaultSettings(
            norm_type, container_name)
        input_method = TemporalRootMeanSquareMethodHelperClass.GetInputMethod(
            container_name)
        output_method = TemporalRootMeanSquareMethodHelperClass.GetOutputMethod(
            container_name)
        container = self.GetContainer(container_name)

        factory = KratosProcessFactory(self.GetModel())
        self.process_list = factory.ConstructListOfProcesses(settings)
        InitializeProcesses(self)

        scalar_list, vec_3d_list, vec_list, mat_list = InitializeContainerArrays(
            container)

        for step in range(0, 12, 2):
            self.model_part.CloneTimeStep(step)
            InitializeModelPartVariables(self.GetModelPart())
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

                analytical_method_scalar = TemporalRootMeanSquareMethodHelperClass.__AnalyticalMethod(
                    norm_type, Kratos.PRESSURE, scalar_list[index])
                analytical_method_vec_3d = TemporalRootMeanSquareMethodHelperClass.__AnalyticalMethod(
                    norm_type, Kratos.VELOCITY, vec_3d_list[index])
                analytical_method_vec = TemporalRootMeanSquareMethodHelperClass.__AnalyticalMethod(
                    norm_type, Kratos.LOAD_MESHES, vec_list[index])
                analytical_method_mat = TemporalRootMeanSquareMethodHelperClass.__AnalyticalMethod(
                    norm_type, Kratos.GREEN_LAGRANGE_STRAIN_TENSOR,
                    mat_list[index])

                if (norm_type == "none"):
                    method_scalar = output_method(item,
                                                  KratosStats.SCALAR_MEAN)
                    method_vec_3d = output_method(item,
                                                  KratosStats.VECTOR_3D_MEAN)
                    method_vec = output_method(item,
                                               Kratos.MATERIAL_PARAMETERS)
                    method_mat = output_method(item,
                                               Kratos.CAUCHY_STRESS_TENSOR)
                else:
                    method_scalar = output_method(item,
                                                  KratosStats.SCALAR_NORM)
                    method_vec_3d = output_method(item,
                                                  KratosStats.VECTOR_3D_NORM)
                    method_vec = output_method(item, Kratos.DENSITY)
                    method_mat = output_method(item, Kratos.VISCOSITY)

                CheckValues(self, analytical_method_scalar, method_scalar, 8)
                CheckValues(self, analytical_method_vec_3d, method_vec_3d, 8)
                CheckValues(self, analytical_method_vec, method_vec, 8)
                CheckValues(self, analytical_method_mat, method_mat, 8)

    @staticmethod
    def __AnalyticalMethod(norm_type, variable, value_array):
        if (norm_type == "none"):
            result = GetInitialVariableValue(variable, "none")
            for item in value_array:
                result += KratosStats.MethodUtilities.RaiseToPower(item, 2) * 2.0
        else:
            result = 0.0
            norm_method = KratosStats.MethodUtilities.GetNormMethod(
                variable, norm_type)
            for item in value_array:
                result += pow(norm_method(item), 2.0) * 2.0

        return KratosStats.MethodUtilities.RaiseToPower(result * (1.0 / max(len(value_array) * 2.0, 1.0)), 0.5)

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
                             "method_name"     : "rootmeansquare",
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
                r'"VECTOR_3D_MEAN", "SCALAR_MEAN", "MATERIAL_PARAMETERS", "CAUCHY_STRESS_TENSOR"'
            )
        else:
            settings_str = settings_str.replace(
                "<OUTPUT_VARIABLES>",
                r'"VECTOR_3D_NORM", "SCALAR_NORM", "DENSITY", "VISCOSITY"')

        return Kratos.Parameters(settings_str)

    @classmethod
    def AddVariables(cls):
        # input variables
        cls.model_part.AddNodalSolutionStepVariable(Kratos.PRESSURE)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.VELOCITY)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.LOAD_MESHES)
        cls.model_part.AddNodalSolutionStepVariable(
            Kratos.GREEN_LAGRANGE_STRAIN_TENSOR)

        # output variables for output_1
        cls.model_part.AddNodalSolutionStepVariable(KratosStats.SCALAR_MEAN)
        cls.model_part.AddNodalSolutionStepVariable(KratosStats.VECTOR_3D_MEAN)
        cls.model_part.AddNodalSolutionStepVariable(
            Kratos.MATERIAL_PARAMETERS)
        cls.model_part.AddNodalSolutionStepVariable(
            Kratos.CAUCHY_STRESS_TENSOR)

        cls.model_part.AddNodalSolutionStepVariable(KratosStats.SCALAR_NORM)
        cls.model_part.AddNodalSolutionStepVariable(KratosStats.VECTOR_3D_NORM)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.DENSITY)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.VISCOSITY)

class TemporalRootMeanSquareMethodTests(
        temporal_statistics_test_case.TemporalStatisticsValueTestCases,
        temporal_statistics_test_case.TemporalStatisticsNormTestCases,
        TemporalRootMeanSquareMethodHelperClass):
    pass


if __name__ == '__main__':
    Kratos.Logger.GetDefaultOutput().SetSeverity(Kratos.Logger.Severity.WARNING)
    import KratosMultiphysics.KratosUnittest as KratosUnittest
    KratosUnittest.main()