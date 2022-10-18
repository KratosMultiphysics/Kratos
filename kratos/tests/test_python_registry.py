import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest

class TestPythonRegistry(KratosUnittest.TestCase):

    def testHasItemNonRegistered(self):
        # Check that provided item is not registered
        self.assertFalse(KratosMultiphysics.Registry.HasItem("FamilyKeyword.NonRegisteredModule.NonRegisteredItem"))

    def testHasItemCpp(self):
        # Check that base Process is registered in its corresponding module as well as in the All block
        self.assertTrue(KratosMultiphysics.Registry.HasItem("Processes.All.Process"))
        self.assertTrue(KratosMultiphysics.Registry.HasItem("Processes.KratosMultiphysics.Process"))

    def testHasItemPython(self):
        # Add a fake entity to the Python registry
        KratosMultiphysics.Registry.AddItem("Processes.KratosMultiphysics.NewProcess", KratosMultiphysics.Process)

        # Check that the fake entity is registered
        self.assertTrue(KratosMultiphysics.Registry.HasItem("Processes.KratosMultiphysics.NewProcess"))

        # Remove the auxiliary testing entity
        KratosMultiphysics.Registry.RemoveItem("Processes.All.NewProcess")
        KratosMultiphysics.Registry.RemoveItem("Processes.KratosMultiphysics.NewProcess")

    def testAddItem(self):
        # Add some fake entities to the Python registry
        KratosMultiphysics.Registry.AddItem("Processes.KratosMultiphysics.NewProcess1", KratosMultiphysics.Process)
        KratosMultiphysics.Registry.AddItem("Processes.KratosMultiphysics.FakeApplication.NewProcess2", KratosMultiphysics.Process)

        # Check that the fake entities are registered
        self.assertTrue(KratosMultiphysics.Registry.HasItem("Processes.All.NewProcess1"))
        self.assertTrue(KratosMultiphysics.Registry.HasItem("Processes.All.NewProcess2"))
        self.assertTrue(KratosMultiphysics.Registry.HasItem("Processes.KratosMultiphysics.NewProcess1"))
        self.assertTrue(KratosMultiphysics.Registry.HasItem("Processes.KratosMultiphysics.FakeApplication.NewProcess2"))

        # Remove the auxiliary testing entities
        KratosMultiphysics.Registry.RemoveItem("Processes.All.NewProcess1")
        KratosMultiphysics.Registry.RemoveItem("Processes.All.NewProcess2")
        KratosMultiphysics.Registry.RemoveItem("Processes.KratosMultiphysics.NewProcess1")
        KratosMultiphysics.Registry.RemoveItem("Processes.KratosMultiphysics.FakeApplication.NewProcess2")

    def testAddItemRepeatedInAllBlock(self):
        # Try to add two entities with same name to different registry modules
        KratosMultiphysics.Registry.AddItem("Processes.KratosMultiphysics.NewProcess", KratosMultiphysics.Process)

        # Check that adding an entity with same family type and name to a different module throws a exception
        with self.assertRaisesRegex(Exception, "Trying to register 'Processes.KratosMultiphysics.FakeApplication.NewProcess' but there is already an item with the same 'NewProcess' name in the 'Processes.All' block."):
            KratosMultiphysics.Registry.AddItem("Processes.KratosMultiphysics.FakeApplication.NewProcess", KratosMultiphysics.Process)

        # Remove the auxiliary testing entity
        KratosMultiphysics.Registry.RemoveItem("Processes.All.NewProcess")
        KratosMultiphysics.Registry.RemoveItem("Processes.KratosMultiphysics.NewProcess")

    def testAddItemAlreadyRegisteredInCpp(self):
        # Try to add an item already present in c++ registry
        with self.assertRaisesRegex(Exception, "Trying to register 'Processes.KratosMultiphysics.Process' but it is already registered."):
            KratosMultiphysics.Registry.AddItem("Processes.KratosMultiphysics.Process", KratosMultiphysics.Process)

    def testAddItemAlreadyRegisteredInPython(self):
        # Add a fake entity to Python registry
        KratosMultiphysics.Registry.AddItem("Processes.KratosMultiphysics.NewProcess", KratosMultiphysics.Process)

        # Try to add an item already present in c++ registry
        with self.assertRaisesRegex(Exception, "Trying to register 'Processes.KratosMultiphysics.NewProcess' but it is already registered."):
            KratosMultiphysics.Registry.AddItem("Processes.KratosMultiphysics.NewProcess", KratosMultiphysics.Process)

        # Remove the auxiliary testing entities
        KratosMultiphysics.Registry.RemoveItem("Processes.All.NewProcess")
        KratosMultiphysics.Registry.RemoveItem("Processes.KratosMultiphysics.NewProcess")

    def testRemoveNonRegistered(self):
        # Check that trying to remove a non-registered item throws an exception
        with self.assertRaisesRegex(Exception, "Asking to remove 'FamilyKeyword.NonRegisteredModule.NonRegisteredItem' non-registered item."):
            KratosMultiphysics.Registry.RemoveItem("FamilyKeyword.NonRegisteredModule.NonRegisteredItem")

    def testRemoveItemCpp(self):
        # Check that trying to remove a c++ registered item throws and exception
        with self.assertRaisesRegex(Exception, "Trying to remove 'Processes.All.Process' from c\+\+ registry. This operation is forbidden."):
            KratosMultiphysics.Registry.RemoveItem("Processes.All.Process")

    def testRemoveItemPython(self):
        # Remove both the module and the all entries one by one
        KratosMultiphysics.Registry.AddItem("FamilyType.Module.ItemName", object())
        KratosMultiphysics.Registry.RemoveItem("FamilyType.All.ItemName")
        KratosMultiphysics.Registry.RemoveItem("FamilyType.Module.ItemName")

        # Remove complete module
        KratosMultiphysics.Registry.AddItem("FamilyType.Module.ItemName1", object())
        KratosMultiphysics.Registry.AddItem("FamilyType.Module.ItemName2", object())
        KratosMultiphysics.Registry.RemoveItem("FamilyType.All")
        KratosMultiphysics.Registry.RemoveItem("FamilyType.Module")

        # Remove all items with same keyword
        # Note that this removes the both the module and the all entries at once
        KratosMultiphysics.Registry.AddItem("FamilyType.Module.ItemName1", object())
        KratosMultiphysics.Registry.AddItem("FamilyType.Module.ItemName2", object())
        KratosMultiphysics.Registry.RemoveItem("FamilyType")

        # Check that registry has no items of FamilyType
        self.assertFalse(KratosMultiphysics.Registry.HasItem("FamilyType"))

    def testGetItemNonRegistered(self):
        # Check that retrieveing a non registered item throws an exception
        with self.assertRaisesRegex(Exception, "Asking to retrieve 'FamilyKeyword.NonRegisteredModule.NonRegisteredItem' non-registered item."):
            KratosMultiphysics.Registry.GetItem("FamilyKeyword.NonRegisteredModule.NonRegisteredItem")

    def testGetItemCppOperation(self):
        # Check the retrieving of a c++ registered operation
        base_operation = KratosMultiphysics.Registry.GetItem("Operations.All.Operation")
        self.assertTrue(isinstance(base_operation, KratosMultiphysics.Operation))

    def testGetItemCppProcess(self):
        # Check the retrieving of a c++ registered process
        base_process = KratosMultiphysics.Registry.GetItem("Processes.All.Process")
        self.assertTrue(isinstance(base_process, KratosMultiphysics.Process))

    def testGetItemPython(self):
        # Add a fake item of the Python registry
        KratosMultiphysics.Registry.AddItem("FamilyType.Module.ItemName", object())

        # Get the test item
        self.assertTrue(isinstance(KratosMultiphysics.Registry.GetItem("FamilyType.All.ItemName"), object))
        self.assertTrue(isinstance(KratosMultiphysics.Registry.GetItem("FamilyType.Module.ItemName"), object))

        # Remove the auxiliary testing entries
        KratosMultiphysics.Registry.RemoveItem("FamilyType.All.ItemName")
        KratosMultiphysics.Registry.RemoveItem("FamilyType.Module.ItemName")

if __name__ == "__main__":
    KratosUnittest.main()
