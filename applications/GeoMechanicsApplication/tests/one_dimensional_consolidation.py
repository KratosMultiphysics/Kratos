import math
import os
import shutil

import KratosMultiphysics as Kratos
import KratosMultiphysics.GeoMechanicsApplication as KratosGeo
import KratosMultiphysics.GeoMechanicsApplication.geomechanics_analysis as analysis

import KratosMultiphysics.KratosUnittest as KratosUnittest

import test_helper
import analytical_solutions


class OneDimensionalConsolidationTestBase(KratosUnittest.TestCase):
    def setUp(self):
        super().setUp()

        self.test_root = test_helper.get_file_path("one_dimensional_consolidation")
        self.test_path = os.path.join(self.test_root, self.get_test_dir_name())

        shutil.rmtree(self.test_path, ignore_errors=True)

        os.makedirs(self.test_path)

        self.number_of_stages = 11
        self.project_parameters_filenames = [f"ProjectParameters_stage{i+1}.json" for i in range(self.number_of_stages)]
        input_filenames = self.project_parameters_filenames + ["MaterialParameters.json", "1D-Consolidationtest.mdpa"]

        for filename in input_filenames:
            shutil.copy(os.path.join(self.test_root, filename), os.path.join(self.test_path, filename))

        self.mid_column_node_ids = [5, 7, 12, 18, 24, 30, 35, 41, 46, 51, 56, 61, 66, 71, 76, 81, 86, 91, 96, 101, 106, 111, 116, 121, 126, 131, 136, 141, 146, 151, 156, 161, 166, 171, 176, 181, 187, 192, 199]
        self.top_node_ids = [197, 198, 199, 200, 201]

        self.end_times = [8640, 17280, 43200, 86400, 172800, 432000, 864000, 1728000, 4320000, 8640000]
        self.t_vs = [0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0]

    def get_test_dir_name(self):
        raise RuntimeError("This base class does not provide a generic test directory name")

class KratosGeoMechanics1DConsolidation(OneDimensionalConsolidationTestBase):
    """
    This class contains benchmark tests which are checked with the analytical solution
    """
    def get_test_dir_name(self):
        return "python"


    def test_1d_consolidation(self):
        """
        test 1D consolidation on elastic soil.

        :return:
        """
        # set stage parameters
        parameters_stages = [None] * self.number_of_stages

        initial_directory = os.getcwd()
        os.chdir(self.test_path)
        for idx, parameter_file_name in enumerate(self.project_parameters_filenames):
            with open(parameter_file_name, 'r') as parameter_file:
                parameters_stages[idx] = Kratos.Parameters(parameter_file.read())

        model = Kratos.Model()
        stages = [analysis.GeoMechanicsAnalysis(model, stage_parameters) for stage_parameters in parameters_stages]

        # run stages and get water pressure/displacement results per stage
        stage_water_pressure = [None] * self.number_of_stages
        stage_displacement   = [None] * self.number_of_stages
        for idx, stage in enumerate(stages):
            stage.Run()
            stage_water_pressure[idx] = test_helper.get_water_pressure(stage)
            displacements = test_helper.get_nodal_variable(stage, KratosGeo.TOTAL_DISPLACEMENT)
            stage_displacement[idx]   = [displacement[1] for displacement in displacements] 

        # get y coords of all the nodes
        post_msh_file_path = os.path.join(self.test_path, "1D-Consolidationtest_stage1.post.msh")
        y_coords = [coord[1] + 1.0 for coord in test_helper.read_coordinates_from_post_msh_file(post_msh_file_path, node_ids=self.mid_column_node_ids)]

        # calculate water pressure analytical solution for all stages and calculate the error
        sample_height = 1.0
        rmse_stages = [None] * (self.number_of_stages - 1)
        for idx, t_v in enumerate(self.t_vs):
            analytical_solution = [analytical_solutions.calculate_relative_water_pressure(y_coord, sample_height, t_v) for y_coord in y_coords]

            # Invert the sign of the water pressures resulting from the numerical solution, to make them match the
            # analytical solution which assumes compressive water pressures to be positive rather than negative
            numerical_solution = [-1.0 * stage_water_pressure[idx + 1][id - 1] for id in self.mid_column_node_ids]

            errors_stage = [rel_p_numerical - rel_p_analytical for rel_p_numerical, rel_p_analytical in
                            zip(numerical_solution, analytical_solution)]
            rmse_stages[idx] = math.sqrt(sum([error * error for error in errors_stage]) / len(errors_stage))

        # assert if average error in all stages is below 1 percent
        accuracy = 0.01
        for idx, rmse_stage in enumerate(rmse_stages):
            self.assertLess(rmse_stage, accuracy, msg=f"RMSE of relative water pressure values in stage {idx+1}")

        # calculate the degree of consolidation analytical solution for all stages and calculate the error
        # Verruijt's notations
        q = -1.0                            # the load
        K = 3.33e2                          # the compression modulus
        G = 5.0e2                           # shear modulus
        m_v = 1/(K + 4/3 * G)               # the compressibility coefficient
        beta = 0.5e-9                       # the compressibility of the water.
        n = 0.3                             # the porosity
        delta_h0 = (-1)*m_v*sample_height*q*n*beta/(m_v+n*beta) # deformation immediately after the application of the load
        delta_h_infinity = (-1)*m_v*sample_height*q # the final deformation

        for idx, t_v in enumerate(self.t_vs):
            rel_displacement = [(-1.0 * stage_displacement[idx + 1][id - 1] - delta_h0) / (delta_h_infinity - delta_h0) for id in
                            self.top_node_ids]
            analytical_degree = analytical_solutions.calculate_degree_of_1d_consolidation(t_v)

            errors_stage = [numerical_degree - analytical_degree for numerical_degree in rel_displacement]
            rmse_stages[idx] = (sum([error ** 2 for error in errors_stage]) / len(errors_stage)) ** 0.5

        # assert if average error in all stages is below 1 percent
        accuracy = 0.01
        for idx, rmse_stage in enumerate(rmse_stages):
            self.assertLess(rmse_stage, accuracy, msg=f"RMSE of degree of consolidation values in stage {idx+1}")

        os.chdir(initial_directory)


class KratosGeoMechanics1DConsolidationCppRoute(OneDimensionalConsolidationTestBase):
    def get_test_dir_name(self):
        return "cpp"


    def test_1d_consolidation(self):
        import KratosMultiphysics.GeoMechanicsApplication.run_geo_settlement as run_geo_settlement

        status = run_geo_settlement.run_stages(self.test_path, self.project_parameters_filenames)
        self.assertEqual(status, 0)

        post_msh_file_path = os.path.join(self.test_path, "1D-Consolidationtest_stage1.post.msh")
        y_coords = [coord[1] + 1.0 for coord in test_helper.read_coordinates_from_post_msh_file(post_msh_file_path, node_ids=self.mid_column_node_ids)]

        sample_height = 1.0
        rmse_stages = []
        for idx, t_v in enumerate(self.t_vs):
            analytical_solution = [analytical_solutions.calculate_relative_water_pressure(y, sample_height, t_v) for y in y_coords]

            output_file_path = os.path.join(self.test_path, f"1D-Consolidationtest_stage{idx+2}.post.res")
            reader = test_helper.GiDOutputFileReader()
            output_data = reader.read_output_from(output_file_path)
            numerical_solution = [-1.0 * pw for pw in reader.nodal_values_at_time("WATER_PRESSURE", self.end_times[idx], output_data, self.mid_column_node_ids)]

            errors_stage = [rel_p_numerical - rel_p_analytical for rel_p_numerical, rel_p_analytical in
                            zip(numerical_solution, analytical_solution)]
            rmse_stages.append(math.sqrt(sum([error * error for error in errors_stage]) / len(errors_stage)))

        # assert if average error in all stages is below 1 percent
        accuracy = 0.01
        for idx, rmse_stage in enumerate(rmse_stages):
            self.assertLess(rmse_stage, accuracy, msg=f"RMSE of relative water pressure values in stage {idx+2}")



if __name__ == '__main__':
    KratosUnittest.main()
