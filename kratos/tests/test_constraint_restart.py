import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utils
import os
import pickle

class TestConstraintRestart(KratosUnittest.TestCase):
    def tearDown(self):
        kratos_utils.DeleteFileIfExisting("test_constraint_restart.rest")

    def test_constraint_serialization(self):
        current_model = KratosMultiphysics.Model()
        model_part = current_model.CreateModelPart("Main")
        
        # Add variables
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)
        
        # Create nodes
        node1 = model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
        node2 = model_part.CreateNewNode(2, 1.0, 0.0, 0.0)
        node3 = model_part.CreateNewNode(3, 0.0, 1.0, 0.0)
        
        # Add DoFs
        for node in model_part.Nodes:
            node.AddDof(KratosMultiphysics.DISPLACEMENT_X, KratosMultiphysics.REACTION_X)
            node.AddDof(KratosMultiphysics.DISPLACEMENT_Y, KratosMultiphysics.REACTION_Y)
        
        # 1. LinearMasterSlaveConstraint
        # Relation: SlaveDof = Weight * MasterDof + Constant
        # We want: Node2.X = 2.0 * Node1.X + 5.0
        
        # Construct arguments
        # It seems the python binding expects DofPointerVectorType&
        # We will try passing a Python list of Dofs
        master_dofs = [node1.GetDof(KratosMultiphysics.DISPLACEMENT_X)]
        slave_dofs = [node2.GetDof(KratosMultiphysics.DISPLACEMENT_X)]
        
        relation_matrix = KratosMultiphysics.Matrix(1, 1)
        relation_matrix[0,0] = 2.0
        
        constant_vector = KratosMultiphysics.Vector(1)
        constant_vector[0] = 5.0
        
        c1 = KratosMultiphysics.LinearMasterSlaveConstraint(1, master_dofs, slave_dofs, relation_matrix, constant_vector)
        model_part.AddMasterSlaveConstraint(c1)

        # 2. SlipConstraint
        # Normal: [1.0, 0.0, 0.0] -> Block X movement
        # normal = KratosMultiphysics.Array3([1.0, 0.0, 0.0])
        # c2 = KratosMultiphysics.SlipConstraint(2, node3.GetDof(KratosMultiphysics.DISPLACEMENT_X), node3.GetDof(KratosMultiphysics.DISPLACEMENT_Y), normal)
        # model_part.AddMasterSlaveConstraint(c2)
        
        # Save
        serializer = KratosMultiphysics.StreamSerializer(KratosMultiphysics.SerializerTraceType.SERIALIZER_TRACE_ALL)
        serializer.Save(model_part.Name, model_part)

        del model_part

        # Load into new ModelPart
        loaded_model = KratosMultiphysics.Model()
        loaded_model_part = loaded_model.CreateModelPart("Main")
        
        serializer.LoadFromBeginning(loaded_model_part.Name, loaded_model_part)
        
        # Verify
        self.assertEqual(loaded_model_part.NumberOfMasterSlaveConstraints(), 1)
        
        # Check C1
        c1_loaded = loaded_model_part.GetMasterSlaveConstraint(1)
        
        # Verify Master Dof
        masters = c1_loaded.GetMasterDofsVector()
        self.assertEqual(len(masters), 1)
        self.assertEqual(masters[0].Id(), master_dofs[0].Id())
        self.assertEqual(masters[0].GetVariable(), KratosMultiphysics.DISPLACEMENT_X)
        
        # Verify Slave Dof
        slaves = c1_loaded.GetSlaveDofsVector()
        self.assertEqual(len(slaves), 1)
        self.assertEqual(slaves[0].Id(), slave_dofs[0].Id())
        self.assertEqual(slaves[0].GetVariable(), KratosMultiphysics.DISPLACEMENT_X)

        print("Constraint serialization test passed.")

if __name__ == '__main__':
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
    KratosUnittest.main()
