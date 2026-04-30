"""Integration tests for PiecewiseLinearMomentCapacityConstitutiveLaw

These tests validate that:
1. The material variables are registered in the GeoMechanics application
2. Material properties can be created and populated
3. The variables are accessible through the Kratos variable system

Note: Detailed moment-curvature response testing is provided in C++ unit tests
(test_piecewise_linear_moment_capacity_constitutive_law.cpp).
"""

import KratosMultiphysics as KM
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.GeoMechanicsApplication import *


class PiecewiseLinearMomentCapacityMaterialTests(KratosUnittest.TestCase):
    """
    Integration tests for the piecewise linear moment capacity constitutive law.
    
    Validates material property system integration with GeoMechanics application.
    """

    def setUp(self):
        """Set up test fixtures"""
        self.model = KM.Model()
        self.model_part = self.model.CreateModelPart("TestModelPart")

    def tearDown(self):
        """Clean up after tests"""
        pass

    def test_material_variables_are_accessible(self):
        """
        Test that the three new material variables are accessible from the GeoMechanicsApplication.
        These are required for the constitutive law to function.
        """
        # Variables should be directly importable from the application
        self.assertIsNotNone(MAX_AXIAL_LOAD_OF_CONSTRUCTION_ELEMENT)
        self.assertIsNotNone(MOMENT_CAPACITY_REDUCTION_FACTOR)
        self.assertIsNotNone(MINIMUM_MOMENT_CAPACITY_FACTOR)

    def test_properties_creation_with_id(self):
        """
        Test that material properties can be created with an ID.
        Properties require an integer ID in Python.
        """
        properties = KM.Properties(0)
        self.assertIsNotNone(properties)
        self.assertEqual(properties.Id, 0)

    def test_properties_can_store_and_retrieve_scalars(self):
        """
        Test that material properties can store and retrieve scalar values.
        """
        properties = KM.Properties(1)
        
        # Set scalar values
        properties.SetValue(MAX_AXIAL_LOAD_OF_CONSTRUCTION_ELEMENT, 10.0)
        properties.SetValue(MOMENT_CAPACITY_REDUCTION_FACTOR, 0.02)
        properties.SetValue(MINIMUM_MOMENT_CAPACITY_FACTOR, 0.30)
        
        # Retrieve and verify
        axial_load = properties.GetValue(MAX_AXIAL_LOAD_OF_CONSTRUCTION_ELEMENT)
        reduction = properties.GetValue(MOMENT_CAPACITY_REDUCTION_FACTOR)
        minimum = properties.GetValue(MINIMUM_MOMENT_CAPACITY_FACTOR)
        
        self.assertAlmostEqual(axial_load, 10.0, places=12)
        self.assertAlmostEqual(reduction, 0.02, places=12)
        self.assertAlmostEqual(minimum, 0.30, places=12)

    def test_properties_scalar_update(self):
        """
        Test that scalar property values can be updated.
        """
        properties = KM.Properties(2)
        
        # Initial value
        properties.SetValue(MAX_AXIAL_LOAD_OF_CONSTRUCTION_ELEMENT, 5.0)
        self.assertAlmostEqual(properties.GetValue(MAX_AXIAL_LOAD_OF_CONSTRUCTION_ELEMENT), 5.0, places=12)
        
        # Update value
        properties.SetValue(MAX_AXIAL_LOAD_OF_CONSTRUCTION_ELEMENT, 15.0)
        self.assertAlmostEqual(properties.GetValue(MAX_AXIAL_LOAD_OF_CONSTRUCTION_ELEMENT), 15.0, places=12)

    def test_properties_in_model_part(self):
        """
        Test that properties can be added to a model part and retrieved.
        This is the typical usage pattern in Kratos simulations.
        """
        # Add properties to model part
        self.model_part.AddProperties(KM.Properties(1))
        
        # Retrieve and configure properties
        properties = self.model_part.GetProperties()[1]
        properties.SetValue(MAX_AXIAL_LOAD_OF_CONSTRUCTION_ELEMENT, 8.0)
        properties.SetValue(MOMENT_CAPACITY_REDUCTION_FACTOR, 0.015)
        
        # Verify
        self.assertAlmostEqual(properties.GetValue(MAX_AXIAL_LOAD_OF_CONSTRUCTION_ELEMENT), 8.0, places=12)
        self.assertAlmostEqual(properties.GetValue(MOMENT_CAPACITY_REDUCTION_FACTOR), 0.015, places=12)

    def test_multiple_property_sets(self):
        """
        Test that multiple property sets can be created and managed independently.
        """
        # Create first property set
        self.model_part.AddProperties(KM.Properties(1))
        props1 = self.model_part.GetProperties()[1]
        props1.SetValue(MAX_AXIAL_LOAD_OF_CONSTRUCTION_ELEMENT, 10.0)
        
        # Create second property set
        self.model_part.AddProperties(KM.Properties(2))
        props2 = self.model_part.GetProperties()[2]
        props2.SetValue(MAX_AXIAL_LOAD_OF_CONSTRUCTION_ELEMENT, 20.0)
        
        # Verify independence
        self.assertAlmostEqual(props1.GetValue(MAX_AXIAL_LOAD_OF_CONSTRUCTION_ELEMENT), 10.0, places=12)
        self.assertAlmostEqual(props2.GetValue(MAX_AXIAL_LOAD_OF_CONSTRUCTION_ELEMENT), 20.0, places=12)


if __name__ == '__main__':
    KratosUnittest.main()
