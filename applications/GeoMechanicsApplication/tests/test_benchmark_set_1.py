# from KratosMultiphysics import * as Kratos

import sys
import os

sys.path.append(os.path.join('..', '..', '..'))
sys.path.append(os.path.join('..', 'python_scripts'))

import KratosMultiphysics as Kratos
import KratosMultiphysics.GeoMechanicsApplication as KratosGeo
import KratosMultiphysics.GeoMechanicsApplication.geomechanics_analysis as analysis

import KratosMultiphysics.KratosUnittest as KratosUnittest

import test_helper

class KratosGeoMechanicsBenchmarkSet1(KratosUnittest.TestCase):
    """
    This class contains benchmark tests which are checked with the analytical solution
    """

    def setUp(self):
        # Code here will be placed BEFORE every test in this TestCase.
        pass

    def tearDown(self):
        # Code here will be placed AFTER every test in this TestCase.
        pass

    def test_benchmark1_1(self):
        """
        In this benchmark Bia-axial shearing test conditions are tested in the Kratos-Geomechanics application.
        In that way, this example can be used to verify elastic deformation with a linear elastic model.
        :return:
        """
        test_name = 'Biaxialshearstresswithlinearelasticmodel'
        file_path = test_helper.get_file_path(os.path.join('.', test_name + '.gid'))

        simulation = test_helper.run_kratos(file_path)

        strain_tensor = test_helper.get_green_lagrange_strain_tensor(simulation)

        # check strain_tensor on all gauss points
        for i in range(len(strain_tensor)):
            for j in range(len(strain_tensor[i])):
                self.assertAlmostEqual(-0.00125, strain_tensor[i][j][0,1])

    # todo
    @KratosUnittest.skip("unit test skipped as it is not ready")
    def test_benchmark1_3(self):
        """
        smooth rigid strip footing on elastic soil
        :return:
        """
        from analytical_solutions import rigid_footing
        import math
        #
        test_name = 'smoothrigidfootingonelasticsoil_2'
        file_path = test_helper.get_file_path(os.path.join('.', test_name + '.gid'))

        simulation = test_helper.run_kratos(file_path)

        cauchy_stress_tensor = test_helper.get_cauchy_stress_tensor(simulation)


        step_size = 0.005
        B = 1
        delta = 0.88
        G = 500
        nu=0.333
        settlement = 0.01
        F = settlement*2*(1+nu)*G/delta

        sigma_v = rigid_footing(0, B, delta, G, nu, settlement)
        print("sigma_v = " + str(sigma_v))

        #
        x = [i*step_size for i in range(int(B/step_size)+1)]
        xmid = [(x[i]+x[i+1])/2 for i in range(len(x)-1)]

       # sigma_v = [2 / math.pi * F / 2 / (B * math.sqrt(1 - (xi / B) ** 2)) for xi in xmid]

        model_part = simulation._list_of_output_processes[0].model_part
        elements = model_part.Elements

        distances = []
        min_distance_elements=[]
        min_distance_idxs = []
        for xi in xmid:
            min_distance = 10**10
            min_distance_element=None
            for element_idx,element in enumerate(elements):

                nodes = element.GetNodes()
                #integration_points = element.GetIntegrationPoints()


                centroid=[0,0,0]
                centre_x = sum([node.X0 for node in nodes])/len(nodes)
                centre_y = sum([node.X0 for node in nodes]) / len(nodes)
                centre_z = sum([node.X0 for node in nodes]) / len(nodes)

                centroid = [[centre_x[i],centre_y[i],centre_z[i]] for i in range(len(centre_x))]
                for node in nodes:
                   centroid[0] =centroid[0]+node.X0
                   centroid[1] = centroid[1] + node.Y0
                   centroid[2] = centroid[2] + node.Z0
                centroid = [x/len(nodes) for x in centroid]
               # for integration_point in integration_points:
                #    distance = compute_distance(integration_point,[xi,0.0,0.0])
                distance = test_helper.compute_distance(centroid, [xi, 0.0, 0.0])
                if distance < min_distance:
                    min_distance = distance
                    min_distance_element = element
                    min_distance_idx = element_idx
            min_distance_elements.append(min_distance_element)
            min_distance_idxs.append(min_distance_idx)
            distances.append(min_distance)

        cauchy_stresses_at_load=[cauchy_stress_tensor[min_distance_idx]for min_distance_idx in min_distance_idxs]
        reaction_force = [(cauchy_stress_at_load[0][3]+cauchy_stress_at_load[1][3]+cauchy_stress_at_load[2][3])/3*step_size for cauchy_stress_at_load in cauchy_stresses_at_load]
        cauchy_stress_tensor[min_distance_idxs]


    # todo
    @KratosUnittest.skip("test should be checked")
    def test_benchmark1_4(self):
        """
        test 1D consolidation on elastic soil.

        :return:
        """
        from analytical_solutions import calculate_1D_consolidation

        # define number of stages
        n_stages = 11

        # get the parameter file names for all stages
        test_name = '1D-Consolidation_all_stages'
        file_path = test_helper.get_file_path(os.path.join('.', test_name))
        parameter_file_names = [os.path.join(file_path, 'ProjectParameters_stage' + str(i + 1) + '.json') for i in
                                range(n_stages)]

        # set stage parameters
        parameters_stages = [None] * n_stages
        os.chdir(file_path)
        for idx, parameter_file_name in enumerate(parameter_file_names):
            with open(parameter_file_name, 'r') as parameter_file:
                parameters_stages[idx] = Kratos.Parameters(parameter_file.read())

        model = Kratos.Model()
        stages = [analysis.GeoMechanicsAnalysis(model, stage_parameters) for stage_parameters in parameters_stages]

        # run stages and get water pressure results per stage
        stage_water_pressure = [None] * n_stages
        for idx, stage in enumerate(stages):
            stage.Run()
            stage_water_pressure[idx] = test_helper.get_water_pressure(stage)

        # get y coords of all the nodes
        coords = test_helper.get_nodal_coordinates(stages[0])
        y_coords = [coord[1] + 1 for coord in coords]

        # calculate analytical solution for all stages and calculate the error
        sample_height = 1
        t_vs = [0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0]
        rmse_stages = [None] * (n_stages-1)
        for idx, t_v in enumerate(t_vs):
            rel_p_stages = [calculate_1D_consolidation(y_coord, sample_height, t_v) * -1 for y_coord in y_coords]
            errors_stage = [stage_water_pressure[idx + 1][node_idx] - rel_p for node_idx, rel_p in
                            enumerate(rel_p_stages)]
            rmse_stages[idx] = (sum([error ** 2 for error in errors_stage]) / len(errors_stage)) ** 0.5

        # assert if average error in all stages is below 1 percent
        accuracy = 0.01
        for rmse_stage in rmse_stages:
            self.assertLess(rmse_stage, accuracy)


    def test_benchmark1_5(self):
        """
        test point load on circular tunnel with beam elements
        :return:
        """
        import math
        from analytical_solutions import calculate_max_deflections_ring, calculate_bending_moments_ring

        # Calculate
        test_name = 'test_tunnel'
        file_path = test_helper.get_file_path(os.path.join('.', test_name + '.gid'))
        simulation = test_helper.run_kratos(file_path)

        # get results
        moments = test_helper.get_moment(simulation)
        displacements = test_helper.get_displacement(simulation)

        max_x_disp = max([abs(displacement[0]) for displacement in displacements])
        max_y_disp = max([abs(displacement[1]) for displacement in displacements])
        z_moments = [moment[2] for moment in moments]
        max_moment, min_moment = max(z_moments), min(z_moments)

        # calculate analytical solution
        point_load = -100
        radius = 0.5
        youngs_modulus = 1e6
        m_inertia = 8.333e-8

        eps_h, eps_v = calculate_max_deflections_ring(point_load, radius, youngs_modulus, m_inertia)

        hor_deform_analytic = abs(eps_h * radius)
        vert_deform_analytic = abs(eps_v * radius * 2)

        max_moment_analytic = calculate_bending_moments_ring(point_load, radius, 0)
        min_moment_analytic = calculate_bending_moments_ring(point_load, radius, math.pi/2)

        # calculate error between analytical solution and numerical solution
        error_max_moment = abs(max_moment_analytic - max_moment)
        error_min_moment = abs(min_moment_analytic - min_moment)

        error_max_x_disp = abs(hor_deform_analytic - max_x_disp)
        error_max_y_disp = abs(vert_deform_analytic - max_y_disp)

        # assert if error is smaller than the precision
        precision = 0.001
        self.assertLess(error_max_moment, abs(max_moment * precision))
        self.assertLess(error_min_moment, abs(min_moment * precision))

        self.assertLess(error_max_x_disp, abs(max_x_disp * precision))
        self.assertLess(error_max_y_disp, abs(max_y_disp * precision))


if __name__ == '__main__':
    suites = KratosUnittest.KratosSuites
    smallSuite = suites['small'] # These tests are executed by the continuous integration tool
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([KratosGeoMechanicsBenchmarkSet1]))
    allSuite = suites['all']
    allSuite.addTests(smallSuite)
    KratosUnittest.runTests(suites)

    #KratosUnittest.runTests(KratosGeoMechanicsBenchmarkSet1())
