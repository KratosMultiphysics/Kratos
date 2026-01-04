
import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest
import numpy

class TestIndexTensorAdaptor(KratosUnittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.current_model = Kratos.Model()
        cls.model_part = cls.current_model.CreateModelPart("Main")
        
        # Create 4 nodes
        n1 = cls.model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
        n2 = cls.model_part.CreateNewNode(2, 1.0, 0.0, 0.0)
        n3 = cls.model_part.CreateNewNode(3, 0.0, 1.0, 0.0)
        n4 = cls.model_part.CreateNewNode(4, 1.0, 1.0, 0.0)
        
        # Create Properties
        props = cls.model_part.CreateNewProperties(1)
        
        # Create 2 Elements (Triangles)
        cls.model_part.CreateNewElement("Element2D3N", 1, [1, 2, 3], props)
        cls.model_part.CreateNewElement("Element2D3N", 2, [2, 4, 3], props)

    def test_IndexTensorAdaptorCollect(self):
        # Test Element Container with IndexTensorAdaptor
        # Assuming IndexTensorAdaptor is exposed under Kratos.TensorAdaptors or similar
        adaptor = Kratos.TensorAdaptors.IndexTensorAdaptor(self.model_part.Elements)
        
        # Check compatibility
        adaptor.Check()
        
        # Collect Data
        adaptor.CollectData()
        data = adaptor.data
        
        expected_indices = numpy.array([
            [1, 2, 3],
            [2, 4, 3]
        ], dtype=numpy.int32)
        
        # Check shape
        self.assertEqual(data.shape[0], 2)
        self.assertEqual(data.shape[1], 3)

        # Check content
        # assertVectorAlmostEqual might not work for 2D arrays directly depending on implementation, 
        # but numpy.array_equal is good for exact int matches.
        self.assertTrue(numpy.array_equal(data, expected_indices), 
                        msg=f"Collected indices mismatch.\nExpected:\n{expected_indices}\nGot:\n{data}")

    def test_IndexTensorAdaptorStoreDataError(self):
        adaptor = Kratos.TensorAdaptors.IndexTensorAdaptor(self.model_part.Elements)
        
        # Test Store Data (Should Throw Error)
        with self.assertRaisesRegex(RuntimeError, "StoreData is not implemented"):
            adaptor.StoreData()

if __name__ == "__main__":
    KratosUnittest.main()
