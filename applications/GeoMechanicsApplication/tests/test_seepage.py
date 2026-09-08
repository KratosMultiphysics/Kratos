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
        This test has a high fixed water pressure at the bottom (higher than a hydrostatic
        profile if the reference coordinate was the top of the column). This forces outflow
        which, due to the seepage boundary, leads to a 0.0 pressure at the seepage boundary
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
        
        # Since the bottom boundary is fixed to a high number, leading to outflow, the 
        # seepage nodes should have pressure = 0
        for node_id, pressure in zip(top_node_ids, water_pressures):
            self.assertAlmostEqual(
                pressure,
                0.0,
            )
        
    def test_three_element_seepage_fixed_bottom_boundary_stop_inflow(self):
        """
        Test with a fixed bottom pressure which is lower than a hydrostatic pressure
        would be when the column is filled. This would induce inflow, but the seepage
        boundary prevents that
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
        nodal_flows = GiDOutputFileReader.nodal_values_at_time(
            "NODAL_WATER_FLOW", 1.0, output_data, top_node_ids
        )
        
        # On the seepage face, the nodal flow should be very small
        for node_id, flow in zip(top_node_ids, nodal_flows):
            self.assertLessEqual(
                abs(flow),
                1.0e-18,
            )
        
    def test_three_element_seepage_flux_bottom_boundary(self):
        """
        Test with forced flux. The nodal outflow is therefore known (the seepage boundary
        allows it) and the pressure should be forces to 0.0 by the seepage boundary
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
        nodal_flows = GiDOutputFileReader.nodal_values_at_time(
            "NODAL_WATER_FLOW", 1.0, output_data, top_node_ids
        )
        # On the seepage face, the nodal flow should be 2.5, due to the forces outflux
        for node_id, flow in zip(top_node_ids, nodal_flows):
            self.assertLessEqual(
                flow,
                2.5,
            )

        water_pressures = GiDOutputFileReader.nodal_values_at_time(
            "WATER_PRESSURE", 1.0, output_data, top_node_ids
        )
        
        # Since the bottom boundary is fixed to a high number, leading to outflow, the 
        # seepage nodes should have pressure = 0
        for node_id, pressure in zip(top_node_ids, water_pressures):
            self.assertAlmostEqual(
                pressure,
                0.0,
            )



if __name__ == "__main__":
    KratosUnittest.main()
