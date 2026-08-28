import os

import KratosMultiphysics.KratosUnittest as KratosUnittest
import test_helper

from KratosMultiphysics.GeoMechanicsApplication.gid_output_file_reader import (
    GiDOutputFileReader,
)


class KratosGeoMechanicsSeepageTests(KratosUnittest.TestCase):
    """
    Test suite for seepage conditions on steady state groundwater flow problems.
    """

    def test_three_element_seepage_fixed_bottom_boundary(self):
        """
        Test case with three vertically stacked 1x1m quadrilateral elements.
        - Top boundary: seepage condition
        - Bottom boundary: normal flux condition (zero flux)
        - Initial condition: uniform zero pressure
        - Material: Van Genuchten retention law
        
        This test verifies that the simulation runs successfully with seepage 
        and normal flux boundary conditions.
        """
        test_name = 'seepage_tests'
        file_path = test_helper.get_file_path(os.path.join('.', test_name, "fixed_bottom_boundary"))
        simulation = test_helper.run_kratos(file_path)
        
        # Read output file
        reader = GiDOutputFileReader()
        output_data = reader.read_output_from(
            os.path.join(file_path, "three_element_seepage_test.post.res")
        )
        
        # Get model part
        model_part = simulation.model.GetModelPart(
            "PorousDomain.porous_computational_model_part"
        )
        
        # Verify that top boundary nodes (y=3.0) have seepage condition applied
        # Nodes 7 and 8 are at y=3.0
        top_node_ids = [7, 8]
        water_pressures = GiDOutputFileReader.nodal_values_at_time(
            "WATER_PRESSURE", 1.0, output_data, top_node_ids
        )
        
        # On a seepage face, pressures should be <= 0 (no positive pressure allowed)
        for node_id, pressure in zip(top_node_ids, water_pressures):
            self.assertLessEqual(
                pressure,
                1.0e-6,
                msg=f"Seepage node {node_id} has positive water pressure {pressure}"
            )
        
        # Verify that simulation converged (implicitly tested by run_kratos)
        # and that we have results for all nodes
        all_node_ids = [node.Id for node in model_part.Nodes]
        self.assertEqual(len(all_node_ids), 8, "Expected 8 nodes in the mesh")
        
        all_pressures = GiDOutputFileReader.nodal_values_at_time(
            "WATER_PRESSURE", 1.0, output_data, all_node_ids
        )
        self.assertEqual(len(all_pressures), 8, "Expected pressure results for all 8 nodes")
        
        # Verify effective saturation is also computed
        all_saturations = GiDOutputFileReader.nodal_values_at_time(
            "EFFECTIVE_SATURATION", 1.0, output_data, all_node_ids
        )
        self.assertEqual(len(all_saturations), 8, "Expected saturation results for all 8 nodes")
        
        # All saturations should be between 0 and 1
        for node_id, saturation in zip(all_node_ids, all_saturations):
            self.assertGreaterEqual(saturation, 0.0, 
                f"Node {node_id} has negative saturation {saturation}")
            self.assertLessEqual(saturation, 1.0, 
                f"Node {node_id} has saturation > 1.0: {saturation}")

    def test_three_element_seepage_fixed_bottom_boundary_stop_inflow(self):
        """
        Test case with three vertically stacked 1x1m quadrilateral elements.
        - Top boundary: seepage condition
        - Bottom boundary: normal flux condition (zero flux)
        - Initial condition: uniform zero pressure
        - Material: Van Genuchten retention law
        
        This test verifies that the simulation runs successfully with seepage 
        and normal flux boundary conditions.
        """
        test_name = 'seepage_tests'
        file_path = test_helper.get_file_path(os.path.join('.', test_name, "fixed_bottom_boundary_stop_inflow"))
        simulation = test_helper.run_kratos(file_path)
        
        # Read output file
        reader = GiDOutputFileReader()
        output_data = reader.read_output_from(
            os.path.join(file_path, "three_element_seepage_test.post.res")
        )
        
        # Get model part
        model_part = simulation.model.GetModelPart(
            "PorousDomain.porous_computational_model_part"
        )
        
        # Verify that top boundary nodes (y=3.0) have seepage condition applied
        # Nodes 7 and 8 are at y=3.0
        top_node_ids = [7, 8]
        water_pressures = GiDOutputFileReader.nodal_values_at_time(
            "WATER_PRESSURE", 1.0, output_data, top_node_ids
        )
        
        # On a seepage face, pressures should be <= 0 (no positive pressure allowed)
        for node_id, pressure in zip(top_node_ids, water_pressures):
            self.assertLessEqual(
                pressure,
                1.0e-6,
                msg=f"Seepage node {node_id} has positive water pressure {pressure}"
            )
        
        # Verify that simulation converged (implicitly tested by run_kratos)
        # and that we have results for all nodes
        all_node_ids = [node.Id for node in model_part.Nodes]
        self.assertEqual(len(all_node_ids), 8, "Expected 8 nodes in the mesh")
        
        all_pressures = GiDOutputFileReader.nodal_values_at_time(
            "WATER_PRESSURE", 1.0, output_data, all_node_ids
        )
        self.assertEqual(len(all_pressures), 8, "Expected pressure results for all 8 nodes")
        
        # Verify effective saturation is also computed
        all_saturations = GiDOutputFileReader.nodal_values_at_time(
            "EFFECTIVE_SATURATION", 1.0, output_data, all_node_ids
        )
        self.assertEqual(len(all_saturations), 8, "Expected saturation results for all 8 nodes")
        
        # All saturations should be between 0 and 1
        for node_id, saturation in zip(all_node_ids, all_saturations):
            self.assertGreaterEqual(saturation, 0.0, 
                f"Node {node_id} has negative saturation {saturation}")
            self.assertLessEqual(saturation, 1.0, 
                f"Node {node_id} has saturation > 1.0: {saturation}")
    
    def test_three_element_seepage_flux_bottom_boundary(self):
        """
        Test case with three vertically stacked 1x1m quadrilateral elements.
        - Top boundary: seepage condition
        - Bottom boundary: normal flux condition (zero flux)
        - Initial condition: uniform zero pressure
        - Material: Van Genuchten retention law
        
        This test verifies that the simulation runs successfully with seepage 
        and normal flux boundary conditions.
        """
        test_name = 'seepage_tests'
        file_path = test_helper.get_file_path(os.path.join('.', test_name, "flux_bottom_boundary"))
        simulation = test_helper.run_kratos(file_path)
        
        # Read output file
        reader = GiDOutputFileReader()
        output_data = reader.read_output_from(
            os.path.join(file_path, "three_element_seepage_test.post.res")
        )
        
        # Get model part
        model_part = simulation.model.GetModelPart(
            "PorousDomain.porous_computational_model_part"
        )
        
        # Verify that top boundary nodes (y=3.0) have seepage condition applied
        # Nodes 7 and 8 are at y=3.0
        top_node_ids = [7, 8]
        water_pressures = GiDOutputFileReader.nodal_values_at_time(
            "WATER_PRESSURE", 1.0, output_data, top_node_ids
        )
        
        # On a seepage face, pressures should be <= 0 (no positive pressure allowed)
        for node_id, pressure in zip(top_node_ids, water_pressures):
            self.assertLessEqual(
                pressure,
                1.0e-6,
                msg=f"Seepage node {node_id} has positive water pressure {pressure}"
            )
        
        # Verify that simulation converged (implicitly tested by run_kratos)
        # and that we have results for all nodes
        all_node_ids = [node.Id for node in model_part.Nodes]
        self.assertEqual(len(all_node_ids), 8, "Expected 8 nodes in the mesh")
        
        all_pressures = GiDOutputFileReader.nodal_values_at_time(
            "WATER_PRESSURE", 1.0, output_data, all_node_ids
        )
        self.assertEqual(len(all_pressures), 8, "Expected pressure results for all 8 nodes")
        
        # Verify effective saturation is also computed
        all_saturations = GiDOutputFileReader.nodal_values_at_time(
            "EFFECTIVE_SATURATION", 1.0, output_data, all_node_ids
        )
        self.assertEqual(len(all_saturations), 8, "Expected saturation results for all 8 nodes")
        
        # All saturations should be between 0 and 1
        for node_id, saturation in zip(all_node_ids, all_saturations):
            self.assertGreaterEqual(saturation, 0.0, 
                f"Node {node_id} has negative saturation {saturation}")
            self.assertLessEqual(saturation, 1.0, 
                f"Node {node_id} has saturation > 1.0: {saturation}")




if __name__ == "__main__":
    KratosUnittest.main()
