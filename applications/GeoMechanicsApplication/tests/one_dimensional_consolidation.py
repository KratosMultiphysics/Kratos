import math
import os

import KratosMultiphysics as Kratos
import KratosMultiphysics.GeoMechanicsApplication as KratosGeo
import KratosMultiphysics.GeoMechanicsApplication.geomechanics_analysis as analysis

import KratosMultiphysics.KratosUnittest as KratosUnittest

import test_helper

class KratosGeoMechanics1DConsolidation(KratosUnittest.TestCase):
    """
    This class contains benchmark tests which are checked with the analytical solution
    """
    def test_1d_consolidation(self):
        """
        test 1D consolidation on elastic soil.

        :return:
        """
        from analytical_solutions import calculate_relative_water_pressure, calculate_degree_of_1d_consolidation
        
        # define number of stages
        n_stages = 11

        # get the parameter file names for all stages
        test_name = 'one_dimensional_consolidation'
        file_path = test_helper.get_file_path(os.path.join('.', test_name))
        parameter_file_names = [os.path.join(file_path, 'ProjectParameters_stage' + str(i + 1) + '.json') for i in
                                range(n_stages)]

        # set stage parameters
        parameters_stages = [None] * n_stages

        initial_directory = os.getcwd()
        os.chdir(file_path)
        for idx, parameter_file_name in enumerate(parameter_file_names):
            with open(parameter_file_name, 'r') as parameter_file:
                parameters_stages[idx] = Kratos.Parameters(parameter_file.read())

        model = Kratos.Model()
        stages = [analysis.GeoMechanicsAnalysis(model, stage_parameters) for stage_parameters in parameters_stages]

        # run stages and get water pressure/displacement results per stage
        stage_water_pressure = [None] * n_stages
        stage_displacement   = [None] * n_stages
        for idx, stage in enumerate(stages):
            stage.Run()
            stage_water_pressure[idx] = test_helper.get_water_pressure(stage)
            displacements = test_helper.get_nodal_variable(stage, KratosGeo.TOTAL_DISPLACEMENT)
            stage_displacement[idx]   = [displacement[1] for displacement in displacements] 

        # get y coords of all the nodes
        node_ids = [5, 7, 12, 18, 24, 30, 35, 41, 46, 51, 56, 61, 66, 71, 76, 81, 86, 91, 96, 101, 106, 111, 116, 121, 126, 131, 136, 141, 146, 151, 156, 161, 166, 171, 176, 181, 187, 192, 199]
        coords = test_helper.get_nodal_coordinates(stages[0])
        y_coords = [1.0 + coords[id - 1][1] for id in node_ids]
        
        t_vs = [0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0]    
        # calculate water pressure analytical solution for all stages and calculate the error
        sample_height = 1.0
        rmse_stages = [None] * (n_stages - 1)
        for idx, t_v in enumerate(t_vs):
            analytical_solution = [calculate_relative_water_pressure(y_coord, sample_height, t_v) for y_coord in y_coords]

            # Compressive water pressures are assumed positive by the analytical solution
            numerical_solution = [-1.0 * stage_water_pressure[idx + 1][id - 1] for id in node_ids]

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

        node_ids = [197, 198, 199, 200, 201]  # nodes at the top edge of the soil column
        for idx, t_v in enumerate(t_vs):
            rel_displacement = [(-1.0 * stage_displacement[idx + 1][id - 1] - delta_h0) / (delta_h_infinity - delta_h0) for id in
                            node_ids]
            analytical_degree_of_consolidation = [calculate_degree_of_1d_consolidation(t_v)
                            for id in node_ids]

            errors_stage = [numerical_degree - analytical_degree for numerical_degree, analytical_degree in zip(rel_displacement, analytical_degree_of_consolidation)]
            rmse_stages[idx] = (sum([error ** 2 for error in errors_stage]) / len(errors_stage)) ** 0.5

        # assert if average error in all stages is below 1 percent
        accuracy = 0.01
        for idx, rmse_stage in enumerate(rmse_stages):
            self.assertLess(rmse_stage, accuracy, msg=f"RMSE of degree of consolidation values in stage {idx+1}")

        os.chdir(initial_directory)

if __name__ == '__main__':
    KratosUnittest.main()
