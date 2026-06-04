from KratosMultiphysics.GeoMechanicsApplication.gid_output_file_reader import GiDOutputFileReader
from KratosMultiphysics import KratosUnittest

class TestCase(KratosUnittest.TestCase):
	def assert_y_displacements_at_time(self, result, node_ids, expected_y, places, time=1.0):
		displacements = GiDOutputFileReader.nodal_values_at_time("DISPLACEMENT", time, result, node_ids=node_ids)
		for displacement in displacements:
			self.assertAlmostEqual(expected_y, displacement[1], places)

	def assert_nodal_values_at_time(self, result, variable_name, expected_values, places, time=1.0):
		for node_id, expected_node_values in enumerate(expected_values, start=1):
			result_values = GiDOutputFileReader.nodal_values_at_time(variable_name, time, result, [node_id])[0]
			self.assertAlmostEqual(result_values, expected_node_values, places)

	def assert_integration_point_tensor_results(self, integration_point_tensors, expected_integration_point_tensor, places, result_name, delta = None):
		for ip_index, ip_tensor in enumerate(integration_point_tensors):
			self.assertVectorAlmostEqual(expected_integration_point_tensor, ip_tensor, places, delta = delta, msg = f"{result_name} components at integration point {ip_index}")

	def assert_integration_point_tensors(self, result, variable_name, expected_tensors, places=None, time=1.0, delta=None):
		for (element_id, ip_index), expected_tensor in expected_tensors.items():
			tensor = GiDOutputFileReader.element_integration_point_values_at_time(variable_name, time, result, [element_id], [ip_index])[0]
			self.assert_integration_point_tensor_results(tensor, expected_tensor, places, variable_name, delta)