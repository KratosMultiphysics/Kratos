import sys
import os

import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.GeoMechanicsApplication.geomechanics_analysis as analysis
import KratosMultiphysics.GeoMechanicsApplication.run_multiple_stages as run_multiple_stages
import test_helper

class KratosGeoMechanicsResetDisplacementTests(KratosUnittest.TestCase):
    """
    This class contains benchmark tests which are checked with the analytical solution
    """
    def test_reset_displacement_truss(self):
        """
        Tests reset displacement in a truss in 4 stages
        stage 1: load is applied
        stage 2: load is kept constant
        stage 3: load is doubled
        stage 4: load is completely removed

        :return:
        """

        # calculate strain
        F = -1e10   # [N]
        E = 2069e8  # [N/m2]
        A = 1       # [m2]

        eps = F/(E*A)

        # get stages
        test_name = 'truss_with_reset_displacement'
        project_path = test_helper.get_file_path(os.path.join('.', test_name))
        n_stages = 4

        cwd = os.getcwd()
        parameter_file_paths = [os.path.join(project_path, f'ProjectParameters_stage{i+1}.json') for i in range(n_stages)]

        # set stage parameters
        parameters_stages = []
        os.chdir(project_path)
        for parameter_file_path in parameter_file_paths:
            with open(parameter_file_path, 'r') as parameter_file:
                parameters_stages.append(KratosMultiphysics.Parameters(parameter_file.read()))

        model = KratosMultiphysics.Model()
        displacement_stages = []
        nodal_coordinates_stages = []

        # run stages and get results
        for stage_parameters in parameters_stages:
            stage = analysis.GeoMechanicsAnalysis(model, stage_parameters)
            stage.Run()
            displacement_stages.append(test_helper.get_displacement(stage))
            nodal_coordinates_stages.append(test_helper. get_nodal_coordinates(stage))

        os.chdir(cwd)

        # Assert
        stage_nr = 1
        for idx,node in enumerate(nodal_coordinates_stages[stage_nr-1]):
            self.assertAlmostEqual(displacement_stages[stage_nr-1][idx][0], eps*node[0], msg = f"u_x at node {idx + 1} in stage {stage_nr}")

        stage_nr = 2
        for idx,node in enumerate(nodal_coordinates_stages[stage_nr-1]):
            self.assertAlmostEqual(displacement_stages[stage_nr-1][idx][0], 0, msg = f"u_x at node {idx + 1} in stage {stage_nr}")

        stage_nr = 3
        for idx,node in enumerate(nodal_coordinates_stages[stage_nr-1]):
            self.assertAlmostEqual(displacement_stages[stage_nr-1][idx][0], eps*node[0], msg = f"u_x at node {idx + 1} in stage {stage_nr}")

        stage_nr = 4
        for idx,node in enumerate(nodal_coordinates_stages[stage_nr-1]):
            self.assertAlmostEqual(displacement_stages[stage_nr-1][idx][0], -2 * eps*node[0], msg = f"u_x at node {idx + 1} in stage {stage_nr}")


    def test_reset_displacement_beam(self):
        """
        Tests reset displacement in a beam in 4 stages
        stage 1: load is applied / reset displacement is true
        stage 2: load is applied / reset displacement is true
        stage 3: load is applied / reset displacement is false
        stage 4: load is removed / reset displacement is false
        """
        project_path = test_helper.get_file_path('geo_beam_with_reset_displacement')
        n_stages = 4
        run_multiple_stages.run_stages(project_path, n_stages)

        # Assert

        # calculate strain
        F = -1e10   # [N]
        E = 2069e8  # [N/m2]
        I = 1       # [m4]
        L = 1       # [m]
        eps = (F*L**3)/(3*E*I)

        reader = test_helper.GiDOutputFileReader()
        output_data = reader.read_output_from(os.path.join(project_path, "geo_beam_with_reset_displacement_stage_1.post.res"))
        time = 1.0
        end_node_id = 11
        y_displacement_at_end_of_beam = reader.nodal_values_at_time("DISPLACEMENT", time, output_data, [end_node_id])[0][1]
        self.assertAlmostEqual(y_displacement_at_end_of_beam, eps * L, places=5)

        output_data = reader.read_output_from(os.path.join(project_path, "geo_beam_with_reset_displacement_stage_2.post.res"))
        time = 2.0
        displacement_vectors = reader.nodal_values_at_time("DISPLACEMENT", time, output_data)
        for u in displacement_vectors:
            self.assertAlmostEqual(u[1], 0.0, places=5)

        output_data = reader.read_output_from(os.path.join(project_path, "geo_beam_with_reset_displacement_stage_3.post.res"))
        time = 3.0
        displacement_vectors = reader.nodal_values_at_time("DISPLACEMENT", time, output_data)
        for u in displacement_vectors:
            self.assertAlmostEqual(u[1], 0.0, places=5)

        output_data = reader.read_output_from(os.path.join(project_path, "geo_beam_with_reset_displacement_stage_4.post.res"))
        time = 4.0
        y_displacement_at_end_of_beam = reader.nodal_values_at_time("DISPLACEMENT", time, output_data, [end_node_id])[0][1]
        self.assertAlmostEqual(y_displacement_at_end_of_beam, -eps * L, places=5)

    def test_reset_displacement_shell_Dirichlet(self):
        """
        Tests reset displacement in a shell, loaded with prescribed Displacement
        Verifies that the prescribed displacement is not erased by reset_displacement.
        load is applied / reset displacement is true, prescribed Z displacement -20 on node 3
        """
        test_name  = 'shell_with_reset_displacement'
        file_path  = test_helper.get_file_path(os.path.join('.', test_name))
        simulation = test_helper.run_kratos(file_path)
        displacement = test_helper.get_displacement(simulation)
        self.assertEqual(displacement[2][2], -20.0)

if __name__ == '__main__':
    KratosUnittest.main()
