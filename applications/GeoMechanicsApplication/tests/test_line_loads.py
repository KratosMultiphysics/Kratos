import os

import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.StructuralMechanicsApplication as KratosStructural
import test_helper
from KratosMultiphysics.project import Project
import importlib
import KratosMultiphysics.GeoMechanicsApplication.run_multiple_stages as run_multiple_stages

class KratosGeoMechanicsLineLoadTests(KratosUnittest.TestCase):
    """
    This class contains benchmark tests to test line loads (in 3D). 
    No analytical solution is defined. This is to check if find neighbour elements work when line elements are applied in 3D.
    """
    def total_reaction_y_from_file(self, output_file_path, time, node_ids):
        output_reader = test_helper.GiDOutputFileReader()
        output_data = output_reader.read_output_from(output_file_path)
        reactions = test_helper.GiDOutputFileReader.nodal_values_at_time("REACTION", time, output_data, node_ids=node_ids)
        return sum([reaction[1] for reaction in reactions])

    def test_line_load_3D2N_hex(self):
        test_name = 'line_load_3D2N_hex'
        parent_name = 'line_load_tests'
        file_path = test_helper.get_file_path(os.path.join(parent_name, test_name + '.gid'))
        simulation = test_helper.run_kratos(file_path)

        displacements = test_helper.get_displacement(simulation)
        y_displacements = [displacement[1] for displacement in displacements]

        node_number = 2
        expected_value = 0.00042943
        self.assertAlmostEqual(expected_value, y_displacements[node_number-1], 6)

    def test_line_load_3D2N_tet(self):
        test_name = 'line_load_3D2N_tet'
        parent_name = 'line_load_tests'
        file_path = test_helper.get_file_path(os.path.join(parent_name, test_name + '.gid'))
        simulation = test_helper.run_kratos(file_path)

        displacements = test_helper.get_displacement(simulation)
        y_displacements = [displacement[1] for displacement in displacements]

        node_number = 4
        expected_value = 0.00034701
        self.assertAlmostEqual(expected_value, y_displacements[node_number-1], 6)

    def test_line_load_3D3N_hex(self):
        test_name = 'line_load_3D3N_hex'
        parent_name = 'line_load_tests'
        file_path = test_helper.get_file_path(os.path.join(parent_name, test_name + '.gid'))
        simulation = test_helper.run_kratos(file_path)

        displacements = test_helper.get_displacement(simulation)
        y_displacements = [displacement[1] for displacement in displacements]

        node_number = 5
        expected_value = 0.0006958
        self.assertAlmostEqual(expected_value, y_displacements[node_number-1], 6)

    def test_line_load_3D3N_tet(self):
        test_name = 'line_load_3D3N_tet'
        parent_name = 'line_load_tests'
        file_path = test_helper.get_file_path(os.path.join(parent_name, test_name + '.gid'))
        simulation = test_helper.run_kratos(file_path)

        displacements = test_helper.get_displacement(simulation)
        y_displacements = [displacement[1] for displacement in displacements]

        node_number = 11
        expected_value = 0.00064506
        self.assertAlmostEqual(expected_value, y_displacements[node_number-1], 6)

    def test_nonuniform_line_load(self):
        test_name = 'non-uniform_line_load'
        parent_name = 'line_load_tests'
        file_path = test_helper.get_file_path(os.path.join(parent_name, test_name))
        simulation = test_helper.run_kratos(file_path)

        bottom_node_ids = [1, 2, 6, 11, 17, 25, 34, 46, 59, 75, 90]
        reactions_y = [reaction[1] for reaction in test_helper.get_nodal_variable(simulation, Kratos.REACTION, bottom_node_ids)]

        self.assertAlmostEqual(20.0, sum(reactions_y), 4)

        line_load_y_by_node_id = [( 91,   0.0),
                                  ( 93,   0.0),
                                  ( 98,  -4.0),
                                  (102,  -9.0),
                                  (106, -10.0),
                                  (111, -10.0),
                                  (114,  -6.0),
                                  (118,  -1.0),
                                  (122,   0.0),
                                  (124,   0.0),
                                  (125,   0.0)]
        top_node_ids = [item[0] for item in line_load_y_by_node_id]
        line_loads_y = [item[1] for item in test_helper.get_nodal_variable(simulation, KratosStructural.LINE_LOAD, top_node_ids)]
        for id_and_load_value, output_load_value in zip(line_load_y_by_node_id, line_loads_y):
            self.assertAlmostEqual(id_and_load_value[1], output_load_value, 4)

    def test_line_loads_in_stages(self):
        test_name = 'line_loads_in_stages'
        parent_name = 'line_load_tests'
        file_path = test_helper.get_file_path(os.path.join(parent_name, test_name, "Legacy"))
        comparison_data = [("test_stage1.post.res", 50.0),
                           ("test_stage2.post.res", 100.0)]
        number_of_stages = len(comparison_data)
        run_multiple_stages.run_stages(file_path, number_of_stages)

        time = 1.0
        bottom_node_ids = [1, 2, 6, 11, 17, 25, 34, 46, 59, 75, 90]
        for output_file_name, expected_total_reaction_y in comparison_data:
            output_file_path = os.path.join(file_path, output_file_name)
            total_reaction_y = self.total_reaction_y_from_file(output_file_path, time, bottom_node_ids)
            self.assertAlmostEqual(total_reaction_y, expected_total_reaction_y, places=3)

    def test_line_loads_with_orchestrator(self):
       test_name = 'line_loads_in_stages'
       parent_name = 'line_load_tests'
       file_path = test_helper.get_file_path(os.path.join(parent_name, test_name, "Orchestrator"))
       comparison_data = [("test_stage1.post.res", 50.0),
                          ("test_stage2.post.res", 100.0)]

       project_parameters_filename = test_helper.get_file_path(os.path.join(file_path, "ProjectParameters.json"))
       # Parse simulation settings and run simulation
       with open(project_parameters_filename, 'r') as parameter_file:
           project_parameters = Kratos.Parameters(parameter_file.read())

       cwd = os.getcwd()
       project = Project(project_parameters)
       os.chdir(file_path)

       orchestrator_reg_entry = Kratos.Registry[project.GetSettings()["orchestrator"]["name"].GetString()]
       orchestrator_module = importlib.import_module(orchestrator_reg_entry["ModuleName"])
       orchestrator_class = getattr(orchestrator_module, orchestrator_reg_entry["ClassName"])
       orchestrator_instance = orchestrator_class(project)
       orchestrator_instance.Run()

       time = 1.0
       bottom_node_ids = [1, 2, 6, 11, 17, 25, 34, 46, 59, 75, 90]
       for output_file_name, expected_total_reaction_y in comparison_data:
           output_file_path = os.path.join(file_path, output_file_name)
           total_reaction_y = self.total_reaction_y_from_file(output_file_path, time, bottom_node_ids)
           self.assertAlmostEqual(total_reaction_y, expected_total_reaction_y, places=3)

       os.chdir(cwd)

if __name__ == '__main__':
    KratosUnittest.main()
