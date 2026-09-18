import KratosMultiphysics as Kratos
import numpy as np
import time

def test_ids_list():
    model = Kratos.Model()
    mp = model.CreateModelPart("test")
    
    # Create some nodes
    num_nodes = 10000
    for i in range(num_nodes):
        mp.CreateNewNode(i + 1, 0.0, 0.0, 0.0)
    
    print(f"Testing IdsList with {num_nodes} nodes...")
    
    # Test PointerVectorSet (Nodes is a PointerVectorSet)
    start_time = time.time()
    ids_list = mp.Nodes.IdsList()
    end_time = time.time()
    
    print(f"IdsList() took {end_time - start_time:.6f} seconds")
    
    # Verification
    expected_ids = np.array([node.Id for node in mp.Nodes], dtype=np.int32)
    
    if np.array_equal(ids_list, expected_ids):
        print("Success: IdsList() matches expected IDs.")
    else:
        print("Error: IdsList() does NOT match expected IDs.")
        print("IdsList:", ids_list[:10])
        print("Expected:", expected_ids[:10])

if __name__ == "__main__":
    test_ids_list()
