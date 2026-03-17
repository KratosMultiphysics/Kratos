import numpy as np
import KratosMultiphysics as KM
#from KratosMultiphysics.TensorAdaptors import GeometriesTensorAdaptor
import KratosMultiphysics.KratosUnittest as KratosUnittest

class TestGeometriesTensorAdaptor(KratosUnittest.TestCase):

    def setUp(self):
        self.current_model = KM.Model()
        self.model_part = self.current_model.CreateModelPart("Main")
        
        # Create Element 1: Standard Triangle
        # Nodes 1(0,0), 2(1,0), 3(0,1)
        self.model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
        self.model_part.CreateNewNode(2, 1.0, 0.0, 0.0)
        self.model_part.CreateNewNode(3, 0.0, 1.0, 0.0)
        
        # Create Element 2: Distorted Triangle
        self.model_part.CreateNewNode(4, 2.0, 2.0, 0.0)
        self.model_part.CreateNewNode(5, 4.0, 2.5, 0.0)
        self.model_part.CreateNewNode(6, 2.5, 4.0, 0.0)
        
        prop = self.model_part.CreateNewProperties(1)
        
        self.elem1 = self.model_part.CreateNewElement("Element2D3N", 1, [1, 2, 3], prop)
        self.elem2 = self.model_part.CreateNewElement("Element2D3N", 2, [4, 5, 6], prop)

    def test_ShapeFunctions_Gauss1(self):
        # Method: Gauss 1 (1 point)
        adaptor = KM.TensorAdaptors.GeometriesTensorAdaptor(self.model_part.Elements,KM.TensorAdaptors.GeometriesTensorAdaptor.DatumType.ShapeFunctions, KM.GeometryData.IntegrationMethod.GI_GAUSS_1)
        adaptor.CollectData()
        
        data = adaptor.data
        # Shape: [2, 1, 3]
        self.assertEqual(data.shape, (2, 1, 3))
        
        # Check values: should be 1/3 for both elements (Partition of Unity)
        for elem_idx in range(2):
            s = 0.0
            for i in range(3):
                val = data[elem_idx, 0, i]
                self.assertAlmostEqual(val, 1.0/3.0)
                s += val
            self.assertAlmostEqual(s, 1.0)
            
        # Specific check for Element 1 (Standard)
        for i in range(3):
            self.assertAlmostEqual(data[0, 0, i], 1.0/3.0)

    def test_ShapeFunctions_Gauss2(self):
        # Method: Gauss 2 (3 points for Triangle)
        adaptor = KM.TensorAdaptors.GeometriesTensorAdaptor(self.model_part.Elements,KM.TensorAdaptors.GeometriesTensorAdaptor.DatumType.ShapeFunctions, KM.GeometryData.IntegrationMethod.GI_GAUSS_2)
        adaptor.CollectData()
        
        data = adaptor.data
        # Shape: [2, 3, 3] -> 3 Gauss points
        self.assertEqual(data.shape, (2, 3, 3))
        
        # Check Partition of Unity for all points
        for elem_idx in range(2):
            for gp in range(3):
                s = sum(data[elem_idx, gp, :])
                self.assertAlmostEqual(s, 1.0)

    def test_ShapeFunctionsDerivatives_Gauss1(self):
        adaptor = KM.TensorAdaptors.GeometriesTensorAdaptor(self.model_part.Elements,KM.TensorAdaptors.GeometriesTensorAdaptor.DatumType.ShapeFunctionDerivatives, KM.GeometryData.IntegrationMethod.GI_GAUSS_1)
        adaptor.CollectData()
        
        data = adaptor.ViewData()
        # Shape: [2, 1, 3, 2]
        self.assertEqual(data.shape, (2, 1, 3, 2))
        
        # Element 1: Standard
        # N1 = 1 - x - y  => dN1/dx = -1, dN1/dy = -1
        self.assertAlmostEqual(data[0, 0, 0, 0], -1.0) # dN1/dx
        self.assertAlmostEqual(data[0, 0, 0, 1], -1.0) # dN1/dy
        
        # Element 2: Distorted. Check Sum(Grad N) = 0
        for elem_idx in range(2):
            sum_dx = sum(data[elem_idx, 0, :, 0])
            sum_dy = sum(data[elem_idx, 0, :, 1])
            self.assertAlmostEqual(sum_dx, 0.0)
            self.assertAlmostEqual(sum_dy, 0.0)

    def test_ShapeFunctionsDerivatives_Gauss2(self):
        adaptor = KM.TensorAdaptors.GeometriesTensorAdaptor(self.model_part.Elements,KM.TensorAdaptors.GeometriesTensorAdaptor.DatumType.ShapeFunctionDerivatives, KM.GeometryData.IntegrationMethod.GI_GAUSS_2)
        adaptor.CollectData()
        
        data = adaptor.ViewData()
        # Shape: [2, 3, 3, 2]
        self.assertEqual(data.shape, (2, 3, 3, 2))
        
        # Verify Partition of Unity Gradients sum to 0
        for elem_idx in range(2):
            for gp in range(3):
                sum_dx = sum(data[elem_idx, gp, :, 0])
                sum_dy = sum(data[elem_idx, gp, :, 1])
                self.assertAlmostEqual(sum_dx, 0.0)
                self.assertAlmostEqual(sum_dy, 0.0)

    def test_Jacobians_Gauss2(self):
        adaptor = KM.TensorAdaptors.GeometriesTensorAdaptor(self.model_part.Elements,KM.TensorAdaptors.GeometriesTensorAdaptor.DatumType.Jacobians, KM.GeometryData.IntegrationMethod.GI_GAUSS_2)
        adaptor.CollectData()
        
        data = adaptor.ViewData()
        # Shape: [2, 3, 2, 2] -> [ielem, igauss, dim, local_dim]
        self.assertEqual(data.shape, (2, 3, 2, 2))
        
        # Element 1 J is constant [[1,0], [0,1]] everywhere?
        # x = (1-L1-L2)*0 + L1*1 + L2*0 = L1
        # y = L2
        # dx/dL1 = 1, dx/dL2 = 0
        # dy/dL1 = 0, dy/dL2 = 1.
        # Yes, constant Jacobian.
        for gp in range(3):
            self.assertAlmostEqual(data[0, gp, 0, 0], 1.0)
            self.assertAlmostEqual(data[0, gp, 0, 1], 0.0)
            self.assertAlmostEqual(data[0, gp, 1, 1], 1.0)
            
        # Element 2 J is constant for linear triangle too.
        # Check that all GPs have same Jacobian for Elem 2.
        for gp in range(1, 3):
             for r in range(2):
                 for c in range(2):
                     self.assertAlmostEqual(data[1, gp, r, c], data[1, 0, r, c])

    def test_Area_Triangle(self):
        # We need Jacobians and Weights to compute Area
        # Area = Sum_gp ( det(J) * W )
        
        # 1. Get Jacobians
        adaptor_J = KM.TensorAdaptors.GeometriesTensorAdaptor(self.model_part.Elements,KM.TensorAdaptors.GeometriesTensorAdaptor.DatumType.Jacobians, KM.GeometryData.IntegrationMethod.GI_GAUSS_2)
        adaptor_J.CollectData()
        J_data = adaptor_J.data # [N_elem, N_gauss, Dim, LocDim]

        # 2. Get Weights
        adaptor_W = KM.TensorAdaptors.GeometriesTensorAdaptor(self.model_part.Elements,KM.TensorAdaptors.GeometriesTensorAdaptor.DatumType.IntegrationWeights, KM.GeometryData.IntegrationMethod.GI_GAUSS_2)
        adaptor_W.CollectData()
        W_data = adaptor_W.data # [N_elem, N_gauss]

        # Calculate Areasº
        # Element 1: (0,0), (1,0), (0,1). Area = 0.5 * Base * Height = 0.5 * 1 * 1 = 0.5
        # Element 2: (2.0, 2.0), (4.0, 2.5), (2.5, 4.0).
        # det(J) = 3.75. Sum(w) = 0.5. Area = 3.75 * 0.5 = 1.875.
        
        calculated_areas = []
        for i_elem in range(2):
            area = 0.0
            n_gauss = W_data.shape[1]
            for i_gauss in range(n_gauss):
                J = J_data[i_elem, i_gauss]
                det_J = np.linalg.det(J)
                
                weight = W_data[i_elem, i_gauss]
                area += det_J * weight
            
            calculated_areas.append(area)
            
        self.assertAlmostEqual(calculated_areas[0], 0.5)
        self.assertAlmostEqual(calculated_areas[1], 1.875)

    def test_Area_CurvedTriangle(self):
        import math
        # Create a new ModelPart for this test to avoid mixing with setUp
        current_model = KM.Model()
        model_part = current_model.CreateModelPart("Curved")
        
        # Nodes as specified:
        # (0,0) (1,0) (0,1) (0.5,0.01) (-0.01,0.5) (sqrt(2)/2,sqrt(2)/2)
        # We must map these to the correct Triangle2D6 order: V1, V2, V3, M12, M23, M31
        # V1=(0,0), V2=(1,0), V3=(0,1)
        # M12 (Edge 1-2): (0.5, 0.01)
        # M23 (Edge 2-3): (sqrt(2)/2, sqrt(2)/2)
        # M31 (Edge 3-1): (-0.01, 0.5)
        
        n1 = model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
        n2 = model_part.CreateNewNode(2, 1.0, 0.0, 0.0)
        n3 = model_part.CreateNewNode(3, 0.0, 1.0, 0.0)
        n4 = model_part.CreateNewNode(4, 0.5, 0.05, 0.0)
        n5 = model_part.CreateNewNode(5, 0.5, 0.5, 0.0)
        n6 = model_part.CreateNewNode(6, -0.05, 0.5, 0.0)
        
        # Create Geometry Triangle2D6
        # Order: V1, V2, V3, M12, M23, M31
        # IDs: 1, 2, 3, 4, 5, 6
        geom = model_part.CreateNewGeometry("Triangle2D6", 1, [1, 2, 3, 4, 5, 6])
        
        integration_methods = [KM.GeometryData.IntegrationMethod.GI_GAUSS_2, KM.GeometryData.IntegrationMethod.GI_GAUSS_3, KM.GeometryData.IntegrationMethod.GI_GAUSS_4]
        # Verify Area = 0.5
        # Use high order integration (Gauss 2) to capture curved boundaries
        for method in integration_methods:        
            # 1. Get Jacobians
            adaptor_J = KM.TensorAdaptors.GeometriesTensorAdaptor(
                model_part.Geometries, 
                KM.TensorAdaptors.GeometriesTensorAdaptor.DatumType.Jacobians, 
                method
            )
            adaptor_J.CollectData()
            J_data = adaptor_J.data
            
            # 2. Get Weights
            adaptor_W = KM.TensorAdaptors.GeometriesTensorAdaptor(
                model_part.Geometries, 
                KM.TensorAdaptors.GeometriesTensorAdaptor.DatumType.IntegrationWeights, 
                method
            )
            adaptor_W.CollectData()
            W_data = adaptor_W.data
            
            area = 0.0
            n_gauss = W_data.shape[1]
            for i_gauss in range(n_gauss):
                J = J_data[0, i_gauss]
                det_J = np.linalg.det(J)
                weight = W_data[0, i_gauss]
                area += det_J * weight
                
            self.assertAlmostEqual(area, 0.5)


if __name__ == '__main__':
    KratosUnittest.main()
