from KratosMultiphysics.GeoMechanicsApplication.gid_output_file_reader import GiDOutputFileReader
from KratosMultiphysics import KratosUnittest

class TestCase(KratosUnittest.TestCase):
	def assert_y_displacements_at_time(self, result, node_ids, expected_y, time, places=None, delta=None):
		displacements = GiDOutputFileReader.nodal_values_at_time("DISPLACEMENT", time, result, node_ids=node_ids)
		for displacement in displacements:
			self.assertAlmostEqual(expected_y, displacement[1], places, delta = delta)

	def assert_nodal_values_at_time(self, result, variable_name, node_ids, expected_values, time, places=None, delta=None):
		for node_id, expected_node_values in zip(node_ids, expected_values):
			result_values = GiDOutputFileReader.nodal_values_at_time(variable_name, time, result, [node_id])[0]
			self.assertAlmostEqual(expected_node_values, result_values, places, delta = delta,
						  	msg = f"There is a difference in the {variable_name} components for node {node_id}")

	def assert_integration_point_tensors(self, result, variable_name, expected_tensors, time, places=None, delta=None):
		for (element_id, ip_index), expected_tensor in expected_tensors.items():
			tensor = GiDOutputFileReader.element_integration_point_values_at_time(variable_name, time, result, [element_id], [ip_index])[0][0]
			self.assertVectorAlmostEqual(expected_tensor, tensor, places, delta = delta, 
							msg = f"There is a difference in the {variable_name} components for element {element_id}, integration point {ip_index}")