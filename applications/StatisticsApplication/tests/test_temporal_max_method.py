import KratosMultiphysics as Kratos
from KratosMultiphysics.process_factory import KratosProcessFactory

import KratosMultiphysics.StatisticsApplication as KratosStats
from KratosMultiphysics.StatisticsApplication.test_utilities import InitializeContainerArrays
from KratosMultiphysics.StatisticsApplication.test_utilities import CheckValues
from KratosMultiphysics.StatisticsApplication.test_utilities import InitializeModelPartVariables
from KratosMultiphysics.StatisticsApplication.test_utilities import InitializeProcesses
from KratosMultiphysics.StatisticsApplication.test_utilities import ExecuteProcessFinalizeSolutionStep

import temporal_statistics_test_case


class TemporalMaxMethodHelperClass(
        temporal_statistics_test_case.TemporalStatisticsTestCase):
    @classmethod
    def AddVariables(cls):
        # input variables
        cls.model_part.AddNodalSolutionStepVariable(Kratos.PRESSURE)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.VELOCITY)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.LOAD_MESHES)
        cls.model_part.AddNodalSolutionStepVariable(
            Kratos.GREEN_LAGRANGE_STRAIN_TENSOR)

        # output variables for output_1
        cls.model_part.AddNodalSolutionStepVariable(KratosStats.SCALAR_NORM)
        cls.model_part.AddNodalSolutionStepVariable(KratosStats.VECTOR_3D_NORM)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.YIELD_STRESS)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.CUTTED_AREA)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.DENSITY)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.NET_INPUT_MATERIAL)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.VISCOSITY)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.WET_VOLUME)

    def RunTemporalStatisticsTest(self, norm_type, container_name):

        settings = TemporalMaxMethodHelperClass.__GetDefaultParameters(
            norm_type, container_name)
        input_method = TemporalMaxMethodHelperClass.GetInputMethod(
            container_name)
        output_method = TemporalMaxMethodHelperClass.GetOutputMethod(
            container_name)
        container = self.GetContainer(container_name)

        factory = KratosProcessFactory(self.GetModel())
        self.process_list = factory.ConstructListOfProcesses(settings)
        InitializeProcesses(self)

        scalar_list, vec_3d_list, vec_list, mat_list = InitializeContainerArrays(
            container)
        step_list = []
        for step in range(0, 12, 2):
            self.model_part.CloneTimeStep(step)
            InitializeModelPartVariables(self.GetModelPart())
            ExecuteProcessFinalizeSolutionStep(self)

            step_list.append(step)
            for index, item in enumerate(container):
                current_scalar = input_method(item, Kratos.PRESSURE)
                current_vector_3d = input_method(item, Kratos.VELOCITY)
                current_vector = input_method(item, Kratos.LOAD_MESHES)
                current_matrix = input_method(
                    item, Kratos.GREEN_LAGRANGE_STRAIN_TENSOR)

                if (step >= 0):
                    scalar_list[index].append(current_scalar)
                    vec_3d_list[index].append(current_vector_3d)
                    vec_list[index].append(current_vector)
                    mat_list[index].append(current_matrix)

                analytical_method_scalar = TemporalMaxMethodHelperClass.__AnalyticalMethod(
                    norm_type, Kratos.PRESSURE, scalar_list[index], step_list)
                analytical_method_vec_3d = TemporalMaxMethodHelperClass.__AnalyticalMethod(
                    norm_type, Kratos.VELOCITY, vec_3d_list[index], step_list)
                analytical_method_vec = TemporalMaxMethodHelperClass.__AnalyticalMethod(
                    norm_type, Kratos.LOAD_MESHES, vec_list[index], step_list)
                analytical_method_mat = TemporalMaxMethodHelperClass.__AnalyticalMethod(
                    norm_type, Kratos.GREEN_LAGRANGE_STRAIN_TENSOR,
                    mat_list[index], step_list)

                mean_method_scalar = output_method(item,
                                                   KratosStats.SCALAR_NORM)
                mean_method_vec_3d = output_method(item,
                                                   KratosStats.VECTOR_3D_NORM)
                mean_method_vec = output_method(item, Kratos.DENSITY)
                mean_method_mat = output_method(item, Kratos.VISCOSITY)
                variance_method_scalar = output_method(item,
                                                       Kratos.YIELD_STRESS)
                variance_method_vec_3d = output_method(item,
                                                       Kratos.CUTTED_AREA)
                variance_method_vec = output_method(item,
                                                    Kratos.NET_INPUT_MATERIAL)
                variance_method_mat = output_method(item, Kratos.WET_VOLUME)

                CheckValues(self, analytical_method_scalar[0],
                            mean_method_scalar, 8)
                CheckValues(self, analytical_method_vec_3d[0],
                            mean_method_vec_3d, 8)
                CheckValues(self, analytical_method_vec[0], mean_method_vec, 8)
                CheckValues(self, analytical_method_mat[0], mean_method_mat, 8)
                CheckValues(self, analytical_method_scalar[1],
                            variance_method_scalar, 8)
                CheckValues(self, analytical_method_vec_3d[1],
                            variance_method_vec_3d, 8)
                CheckValues(self, analytical_method_vec[1],
                            variance_method_vec, 8)
                CheckValues(self, analytical_method_mat[1],
                            variance_method_mat, 8)

    @staticmethod
    def __AnalyticalMethod(norm_type, variable, value_array, step_list):
        result_min = -1e+30
        result_min_time = 0.0
        norm_method = KratosStats.MethodUtilities.GetNormMethod(
            variable, norm_type)
        for index, item in enumerate(value_array):
            value = norm_method(item)
            if (result_min < value):
                result_min = value
                result_min_time = step_list[index]

        return result_min, result_min_time

    @staticmethod
    def __GetDefaultParameters(norm_type, container_name):
        settings_str = r'''
        [
            {
                "kratos_module" : "KratosMultiphysics.StatisticsApplication",
                "python_module" : "temporal_statistics_process",
                "Parameters" : {
                    "model_part_name"                : "test_model_part",
                    "input_variable_settings" : [
                        {
                             "method_name"     : "max",
                             "norm_type"       : "<TEST_NORM_TYPE>",
                             "container"       : "<TEST_CONTAINER>",
                             "echo_level"      : 0,
                             "method_settings" : {
                                "input_variables"  : ["VELOCITY", "PRESSURE", "LOAD_MESHES", "GREEN_LAGRANGE_STRAIN_TENSOR"],
                                "output_variables"     : ["VECTOR_3D_NORM", "SCALAR_NORM", "DENSITY", "VISCOSITY"],
                                "output_time_step_variables" : ["CUTTED_AREA", "YIELD_STRESS", "NET_INPUT_MATERIAL", "WET_VOLUME"]
                             }
                        }
                    ],
                    "statistics_start_point_control_variable_name" : "TIME",
                    "statistics_start_point_control_value"         : 0.0
                }
            }
        ]'''
        settings_str = settings_str.replace("<TEST_NORM_TYPE>", norm_type)
        settings_str = settings_str.replace("<TEST_CONTAINER>", container_name)
        return Kratos.Parameters(settings_str)


class TemporalMaxMethodTests(
        temporal_statistics_test_case.TemporalStatisticsNormTestCases,
        TemporalMaxMethodHelperClass):
    pass


if __name__ == '__main__':
    Kratos.Logger.GetDefaultOutput().SetSeverity(Kratos.Logger.Severity.WARNING)
    import KratosMultiphysics.KratosUnittest as KratosUnittest
    KratosUnittest.main()