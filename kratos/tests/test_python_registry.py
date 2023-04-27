import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.python_registry import RegistryContext

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
        KratosMultiphysics.Registry.AddItem("Processes.KratosMultiphysics.NewProcess", KratosMultiphysics.Process())

        # Check that the fake entity is registered
        self.assertTrue(KratosMultiphysics.Registry.HasItem("Processes.KratosMultiphysics.NewProcess"))

        # Remove the auxiliary testing entity from the Python registry
        KratosMultiphysics.Registry.RemoveItem("Processes")

    def testNumberOfItems(self):
        # Add some fake entities to the Python registr
        KratosMultiphysics.Registry.AddItem("FakeEntities.MyModule1.MyObject11", object())
        KratosMultiphysics.Registry.AddItem("FakeEntities.MyModule1.MyObject12", object())
        KratosMultiphysics.Registry.AddItem("FakeEntities.MyModule2.MyObject21", object())

        # Check the numer of items
        self.assertEqual(KratosMultiphysics.Registry.NumberOfItems("FakeEntities.All"), 3)
        self.assertEqual(KratosMultiphysics.Registry.NumberOfItems("FakeEntities.MyModule1"), 2)
        self.assertEqual(KratosMultiphysics.Registry.NumberOfItems("FakeEntities.MyModule2"), 1)
        self.assertEqual(KratosMultiphysics.Registry.NumberOfItems("FakeEntities.All", RegistryContext.CPP), 0)
        self.assertEqual(KratosMultiphysics.Registry.NumberOfItems("FakeEntities.All", RegistryContext.PYTHON), 3)

        # Remove the auxiliary testing tentities from the Python registry
        KratosMultiphysics.Registry.RemoveItem("FakeEntities")


    def testAddItem(self):
        # Add some fake entities to the Python registry
        KratosMultiphysics.Registry.AddItem("Processes.KratosMultiphysics.NewProcess1", KratosMultiphysics.Process())
        KratosMultiphysics.Registry.AddItem("Processes.KratosMultiphysics.FakeApplication.NewProcess2", KratosMultiphysics.Process())

        # Check that the fake entities are registered
        self.assertTrue(KratosMultiphysics.Registry.HasItem("Processes.All.NewProcess1"))
        self.assertTrue(KratosMultiphysics.Registry.HasItem("Processes.All.NewProcess2"))
        self.assertTrue(KratosMultiphysics.Registry.HasItem("Processes.KratosMultiphysics.NewProcess1"))
        self.assertTrue(KratosMultiphysics.Registry.HasItem("Processes.KratosMultiphysics.FakeApplication.NewProcess2"))

        # Remove the auxiliary testing entities from the Python registry
        KratosMultiphysics.Registry.RemoveItem("Processes")

    def testAddItemRepeatedInAllBlock(self):
        # Try to add two entities with same name to different registry modules
        KratosMultiphysics.Registry.AddItem("Processes.KratosMultiphysics.NewProcess", KratosMultiphysics.Process())

        # Check that adding an entity with same family type and name to a different module throws a exception
        with self.assertRaisesRegex(Exception, "Trying to register 'Processes.KratosMultiphysics.FakeApplication.NewProcess' but there is already an item with the same 'NewProcess' name in the 'Processes.All' block."):
            KratosMultiphysics.Registry.AddItem("Processes.KratosMultiphysics.FakeApplication.NewProcess", KratosMultiphysics.Process())

        # Remove the auxiliary testing entity from the Python registry
        KratosMultiphysics.Registry.RemoveItem("Processes")

    def testAddItemAlreadyRegisteredInCpp(self):
        # Try to add an item already present in the c++ registry
        with self.assertRaisesRegex(Exception, "Trying to register 'Processes.KratosMultiphysics.Process' but it is already registered."):
            KratosMultiphysics.Registry.AddItem("Processes.KratosMultiphysics.Process", KratosMultiphysics.Process())

    def testAddItemAlreadyRegisteredInPython(self):
        # Add a fake entity to Python registry
        KratosMultiphysics.Registry.AddItem("Processes.KratosMultiphysics.NewProcess", KratosMultiphysics.Process())

        # Try to add an item already present in the Python registry
        with self.assertRaisesRegex(Exception, "Trying to register 'Processes.KratosMultiphysics.NewProcess' but it is already registered."):
            KratosMultiphysics.Registry.AddItem("Processes.KratosMultiphysics.NewProcess", KratosMultiphysics.Process())

        # Remove the auxiliary testing entities from the Python registry
        KratosMultiphysics.Registry.RemoveItem("Processes")

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
            KratosMultiphysics.Registry["FamilyKeyword.NonRegisteredModule.NonRegisteredItem"]

    def testGetItemCppOperation(self):
        # Check the retrieving of a c++ registered operation
        base_operation = KratosMultiphysics.Registry["Operations.All.Operation.Prototype"]
        self.assertTrue(isinstance(base_operation, KratosMultiphysics.Operation))

    def testGetItemCppProcess(self):
        # Check the retrieving of a c++ registered process
        base_process = KratosMultiphysics.Registry["Processes.All.Process.Prototype"]
        self.assertTrue(isinstance(base_process, KratosMultiphysics.Process))

    def testGetItemPython(self):
        # Add a fake item of the Python registry
        KratosMultiphysics.Registry.AddItem("FamilyType.Module.ItemName", object())

        # Get the test item
        self.assertTrue(isinstance(KratosMultiphysics.Registry["FamilyType.All.ItemName"], object))
        self.assertTrue(isinstance(KratosMultiphysics.Registry["FamilyType.Module.ItemName"], object))

        # Remove the auxiliary testing entries from the Python registry
        KratosMultiphysics.Registry.RemoveItem("FamilyType")

    def testKeys(self):
        # Add a fake items to the Python registry
        KratosMultiphysics.Registry.AddItem("Processes.Module.PythonProcessInModule", object())
        KratosMultiphysics.Registry.AddItem("Processes.KratosMultiphysics.PythonProcessInKratosMultiphysics", object())

        # Check current Processes keys
        expected_proc_keys_cpp = ["All","KratosMultiphysics"]
        expected_proc_keys_py = ["All","Module","KratosMultiphysics"]
        proc_keys_py, proc_keys_cpp = KratosMultiphysics.Registry.keys("Processes")
        for proc_key in proc_keys_py : self.assertTrue(proc_key in expected_proc_keys_py)
        for proc_key in proc_keys_cpp : self.assertTrue(proc_key in expected_proc_keys_cpp)

        # Check current Processes.KratosMultiphysics keys
        expected_proc_keys_cpp = ["Process","OutputProcess"]
        expected_proc_keys_py = ["PythonProcessInKratosMultiphysics"]
        proc_keys_py, proc_keys_cpp = KratosMultiphysics.Registry.keys("Processes.KratosMultiphysics")
        for proc_key in proc_keys_py : self.assertTrue(proc_key in expected_proc_keys_py)
        for proc_key in proc_keys_cpp : self.assertTrue(proc_key in expected_proc_keys_cpp)

        # Check current Processes.All keys
        expected_proc_keys_cpp = ["Process","OutputProcess"]
        expected_proc_keys_py = ["PythonProcessInKratosMultiphysics","PythonProcessInModule"]
        proc_keys_py, proc_keys_cpp = KratosMultiphysics.Registry.keys("Processes.All")
        for proc_key in proc_keys_py : self.assertTrue(proc_key in expected_proc_keys_py)
        for proc_key in proc_keys_cpp : self.assertTrue(proc_key in expected_proc_keys_cpp)

        # Check that trying to retrieve the keys of a value items throws a exception
        with self.assertRaisesRegex(Exception, "Asking for the keys of 'Processes.KratosMultiphysics.Process.Prototype'. 'Processes.KratosMultiphysics.Process.Prototype' item has no subitems."):
            KratosMultiphysics.Registry.keys("Processes.KratosMultiphysics.Process.Prototype")

        # Remove the auxiliary testing entries from the Python registry
        KratosMultiphysics.Registry.RemoveItem("Processes")

    # #FIXME: This way of checking the iteration will most probably crash once we add more stuff to the registry
    # #FIXME: Most probably we should check that we're iterating more than one item or something of this sort
    # def testIteration(self):
    #     KratosMultiphysics.Registry.AddItem("Processes.KratosMultiphysics.KratosApplication.PythonProcess", KratosMultiphysics.Process())
    #     KratosMultiphysics.Registry.AddItem("PythonRootItem.PythonSubItem.PythonSubSubItem", object())

    #     root_items_keys = ["Operations","Processes","PythonRootItem"]
    #     sub_items_keys = ["All","KratosMultiphysics","PythonSubItem"]
    #     sub_sub_items_keys = ["PythonSubSubItem","KratosApplication","PythonProcess","Operation","Process"]
    #     sub_sub_sub_items_keys = ["PythonProcess","ModulePath","Prototype","CreateFunction"]
    #     for item in KratosMultiphysics.Registry:
    #         self.assertTrue(item in root_items_keys)
    #         if not KratosMultiphysics.Registry.HasValue(item):
    #             for subitem in KratosMultiphysics.Registry[item]:
    #                 self.assertTrue(subitem in sub_items_keys)
    #                 if not KratosMultiphysics.Registry.HasValue(f"{item}.{subitem}"):
    #                     for subsubitem in KratosMultiphysics.Registry[f"{item}.{subitem}"]:
    #                         self.assertTrue(subsubitem in sub_sub_items_keys)
    #                     if not KratosMultiphysics.Registry.HasValue(f"{item}.{subitem}.{subsubitem}"):
    #                         for subsubsubitem in KratosMultiphysics.Registry[f"{item}.{subitem}.{subsubitem}"]:
    #                             self.assertTrue(subsubsubitem in sub_sub_sub_items_keys)

    #     # Remove the auxiliary testing entries from the Python registry
    #     KratosMultiphysics.Registry.RemoveItem("Processes")
    #     KratosMultiphysics.Registry.RemoveItem("PythonRootItem")

    def testDecorator(self):
        # Auxiliary process class to be used in the testing
        @KratosMultiphysics.RegisterPrototype("Processes.KratosMultiphysics")
        class FooProcess(KratosMultiphysics.Process):
            def __init__(self, a):
                super().__init__()
                self.a = a

            def getA(self):
                return self.a

        # Assert that the decorator-based registry works
        self.assertTrue(KratosMultiphysics.Registry.HasItem("Processes.All.FooProcess"))
        self.assertTrue(KratosMultiphysics.Registry.HasItem("Processes.KratosMultiphysics.FooProcess"))
        self.assertEqual(KratosMultiphysics.Registry["Processes.KratosMultiphysics.FooProcess.Prototype"](10).getA(), 10)

        # Remove the testing entry from the Python registry
        KratosMultiphysics.Registry.RemoveItem("Processes")

if __name__ == "__main__":
    KratosUnittest.main()
