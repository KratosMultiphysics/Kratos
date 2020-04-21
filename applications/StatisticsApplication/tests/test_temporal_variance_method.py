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


class TemporalVarianceMethodHelperClass(
        temporal_statistics_test_case.TemporalStatisticsTestCase):
    def RunTemporalStatisticsTest(self, norm_type, container_name):

        settings = TemporalVarianceMethodHelperClass.__GetDefaultSettings(
            norm_type, container_name)
        input_method = TemporalVarianceMethodHelperClass.GetInputMethod(
            container_name)
        output_method = TemporalVarianceMethodHelperClass.GetOutputMethod(
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

                analytical_method_scalar = TemporalVarianceMethodHelperClass.__AnalyticalMethod(
                    norm_type, Kratos.PRESSURE, scalar_list[index])
                analytical_method_vec_3d = TemporalVarianceMethodHelperClass.__AnalyticalMethod(
                    norm_type, Kratos.VELOCITY, vec_3d_list[index])
                analytical_method_vec = TemporalVarianceMethodHelperClass.__AnalyticalMethod(
                    norm_type, Kratos.LOAD_MESHES, vec_list[index])
                analytical_method_mat = TemporalVarianceMethodHelperClass.__AnalyticalMethod(
                    norm_type, Kratos.GREEN_LAGRANGE_STRAIN_TENSOR,
                    mat_list[index])

                if (norm_type == "none"):
                    mean_method_scalar = output_method(item,
                                                       KratosStats.SCALAR_MEAN)
                    mean_method_vec_3d = output_method(
                        item, KratosStats.VECTOR_3D_MEAN)
                    mean_method_vec = output_method(item,
                                                    Kratos.MATERIAL_PARAMETERS)
                    mean_method_mat = output_method(
                        item, Kratos.CAUCHY_STRESS_TENSOR)
                    variance_method_scalar = output_method(
                        item, KratosStats.SCALAR_VARIANCE)
                    variance_method_vec_3d = output_method(
                        item, KratosStats.VECTOR_3D_VARIANCE)
                    variance_method_vec = output_method(
                        item, Kratos.ELEMENTAL_DISTANCES)
                    variance_method_mat = output_method(
                        item, Kratos.LOCAL_INERTIA_TENSOR)
                else:
                    mean_method_scalar = output_method(item,
                                                       KratosStats.SCALAR_NORM)
                    mean_method_vec_3d = output_method(
                        item, KratosStats.VECTOR_3D_NORM)
                    mean_method_vec = output_method(item, Kratos.DENSITY)
                    mean_method_mat = output_method(item, Kratos.VISCOSITY)
                    variance_method_scalar = output_method(
                        item, Kratos.YIELD_STRESS)
                    variance_method_vec_3d = output_method(
                        item, Kratos.CUTTED_AREA)
                    variance_method_vec = output_method(
                        item, Kratos.NET_INPUT_MATERIAL)
                    variance_method_mat = output_method(
                        item, Kratos.WET_VOLUME)

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
    def __AnalyticalMethod(norm_type, variable, value_array):
        if (norm_type == "none"):
            result_mean = GetInitialVariableValue(variable, "none")
            result_variance = GetInitialVariableValue(variable, "none")
            for item in value_array:
                result_mean += item * 2.0
                result_variance += KratosStats.MethodUtilities.RaiseToPower(
                    item, 2) * 2.0
        else:
            result_mean = 0.0
            result_variance = 0.0
            norm_method = KratosStats.MethodUtilities.GetNormMethod(
                variable, norm_type)
            for item in value_array:
                value = norm_method(item)
                result_mean += value * 2.0
                result_variance += pow(value, 2) * 2.0

        result_mean = result_mean * (1.0 / max(len(value_array) * 2.0, 1.0))
        result_variance = result_variance * (
            1.0 / max(len(value_array) * 2.0, 1.0)
        ) - KratosStats.MethodUtilities.RaiseToPower(result_mean, 2)
        return result_mean, result_variance

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
                             "method_name"     : "variance",
                             "norm_type"       : "<TEST_NORM_TYPE>",
                             "container"       : "<TEST_CONTAINER>",
                             "echo_level"      : 0,
                             "method_settings" : {
                                "input_variables"  : ["VELOCITY", "PRESSURE", "LOAD_MESHES", "GREEN_LAGRANGE_STRAIN_TENSOR"],
                                "output_mean_variables"     : [<OUTPUT_VARIABLES_1>],
                                "output_variance_variables" : [<OUTPUT_VARIABLES_2>]
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
                "<OUTPUT_VARIABLES_1>",
                r'"VECTOR_3D_MEAN", "SCALAR_MEAN", "MATERIAL_PARAMETERS", "CAUCHY_STRESS_TENSOR"'
            )
            settings_str = settings_str.replace(
                "<OUTPUT_VARIABLES_2>",
                r'"VECTOR_3D_VARIANCE", "SCALAR_VARIANCE", "ELEMENTAL_DISTANCES", "LOCAL_INERTIA_TENSOR"'
            )
        else:
            settings_str = settings_str.replace(
                "<OUTPUT_VARIABLES_1>",
                r'"VECTOR_3D_NORM", "SCALAR_NORM", "DENSITY", "VISCOSITY"')
            settings_str = settings_str.replace(
                "<OUTPUT_VARIABLES_2>",
                r'"CUTTED_AREA", "YIELD_STRESS", "NET_INPUT_MATERIAL", "WET_VOLUME"'
            )

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
            KratosStats.SCALAR_VARIANCE)
        cls.model_part.AddNodalSolutionStepVariable(
            KratosStats.VECTOR_3D_VARIANCE)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.MATERIAL_PARAMETERS)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.ELEMENTAL_DISTANCES)
        cls.model_part.AddNodalSolutionStepVariable(
            Kratos.CAUCHY_STRESS_TENSOR)
        cls.model_part.AddNodalSolutionStepVariable(
            Kratos.LOCAL_INERTIA_TENSOR)

        cls.model_part.AddNodalSolutionStepVariable(KratosStats.SCALAR_NORM)
        cls.model_part.AddNodalSolutionStepVariable(KratosStats.VECTOR_3D_NORM)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.YIELD_STRESS)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.CUTTED_AREA)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.DENSITY)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.NET_INPUT_MATERIAL)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.VISCOSITY)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.WET_VOLUME)


class TemporalVarianceMethodTests(
        temporal_statistics_test_case.TemporalStatisticsValueTestCases,
        temporal_statistics_test_case.TemporalStatisticsNormTestCases,
        TemporalVarianceMethodHelperClass):
    pass


if __name__ == '__main__':
    Kratos.Logger.GetDefaultOutput().SetSeverity(Kratos.Logger.Severity.WARNING)
    import KratosMultiphysics.KratosUnittest as KratosUnittest
    KratosUnittest.main()