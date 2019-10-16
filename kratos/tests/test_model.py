from __future__ import print_function, absolute_import, division

import sys

import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utils

# Use cPickle on Python 2.7 (Note that only the cPickle module is supported on Python 2.7)
# Source: https://pybind11.readthedocs.io/en/stable/advanced/classes.html
pickle_message = ""
try:
    import cPickle as pickle
    have_pickle_module = True
except ImportError:
    if sys.version_info > (3, 0):
        try:
            import pickle
            have_pickle_module = True
        except ImportError:
            have_pickle_module = False
            pickle_message = "No pickle module found"
    else:
        have_pickle_module = False
        pickle_message = "No valid pickle module found"

class TestModel(KratosUnittest.TestCase):

    def test_model(self):
        current_model = KratosMultiphysics.Model()

        model_part = current_model.CreateModelPart("Main")
        model_part.CreateSubModelPart("Inlets")
        model_part.CreateSubModelPart("Temp")
        outlet = model_part.CreateSubModelPart("Outlet")

        aaa = current_model["Main.Outlet"].CreateSubModelPart("aaa")

        self.assertEqual(aaa, current_model["aaa"]) #search by flat name - should be eventually deprecated

        #check that a meaningful error is thrown
        with self.assertRaisesRegex(RuntimeError, "Error: The ModelPart named : \"abc\" was not found either as root-ModelPart or as a flat name. The total input string was \"abc\""):
            current_model["abc"]

        #check that a meaningful error is thrown
        with self.assertRaisesRegex(RuntimeError, "The ModelPart named : \"aaa\" was not found as SubModelPart of : \"Inlets\". The total input string was \"Main.Inlets.aaa\""):
            current_model["Main.Inlets.aaa"]

        #here i create a model part in the lowest level and i check that the other model parts have it
        current_model["Main.Outlet.aaa"].CreateNewNode(1,0.0,0.0,0.0)
        self.assertEqual(len(model_part.Nodes), 1)
        self.assertEqual(len(outlet.Nodes), 1)
        self.assertEqual(len(aaa.Nodes), 1)

        #self.assertEqual(current_model["Main"], model_part )
        #self.assertEqual(current_model["Main.Outlet"].Info(), outlet.Info() )
        #self.assertEqual(current_model["Main.Outlet.aaa"].Info(), aaa.Info() )

        self.assertTrue(current_model.HasModelPart("Main"))
        self.assertTrue(current_model.HasModelPart("Main.Outlet"))
        self.assertFalse(current_model.HasModelPart("Outlet"))

    def _create_and_save_model(self,file_name, serializer_flag):
        current_model = KratosMultiphysics.Model()

        model_part = current_model.CreateModelPart("Main")
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.TEMPERATURE)
        model_part.CreateSubModelPart("Inlets")
        model_part.CreateSubModelPart("Temp")
        model_part.CreateNewNode(1,0.0,0.0,0.0)
        other = current_model.CreateModelPart("Other")
        other.AddNodalSolutionStepVariable(KratosMultiphysics.PRESSURE)
        other.CreateNewNode(1,0.0,0.0,0.0)

        KratosMultiphysics.FileSerializer(file_name, serializer_flag).Save("ModelSerialization",current_model)

    def test_model_serialization(self):

        file_name = "model_serialization"
        serializer_flag = KratosMultiphysics.SerializerTraceType.SERIALIZER_NO_TRACE

        self._create_and_save_model(file_name, serializer_flag)

        loaded_model = KratosMultiphysics.Model()
        KratosMultiphysics.FileSerializer(file_name, serializer_flag).Load("ModelSerialization",loaded_model)

        self.assertTrue(loaded_model["Main"].HasNodalSolutionStepVariable(KratosMultiphysics.TEMPERATURE))
        self.assertTrue(loaded_model["Other"].HasNodalSolutionStepVariable(KratosMultiphysics.PRESSURE))

        self.assertTrue(loaded_model.HasModelPart("Main.Inlets"))
        self.assertTrue(loaded_model.HasModelPart("Main.Temp"))
        self.assertTrue(1 in loaded_model["Main"].Nodes)
        self.assertTrue(1 in loaded_model["Other"].Nodes)

        kratos_utils.DeleteFileIfExisting(file_name + ".rest")

    @KratosUnittest.skipUnless(have_pickle_module, "Pickle module error: " + pickle_message)
    def test_model_serialization_with_pickling(self):
        current_model = KratosMultiphysics.Model()

        model_part = current_model.CreateModelPart("Main")
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.TEMPERATURE)
        model_part.CreateSubModelPart("Inlets")
        model_part.CreateSubModelPart("Temp")
        model_part.CreateNewNode(1,0.0,0.0,0.0)
        other = current_model.CreateModelPart("Other")
        other.AddNodalSolutionStepVariable(KratosMultiphysics.PRESSURE)
        other.CreateNewNode(1,0.0,0.0,0.0)

        serializer = KratosMultiphysics.StreamSerializer()
        serializer.Save("ModelSerialization",current_model)
        del(current_model)

        #pickle dataserialized_data
        pickled_data = pickle.dumps(serializer, protocol=2) # Second argument is the protocol and is NECESSARY (according to pybind11 docs)
        del(serializer)

        #overwrite the old serializer with the unpickled one
        serializer = pickle.loads(pickled_data)

        loaded_model = KratosMultiphysics.Model()
        serializer.Load("ModelSerialization",loaded_model)

        self.assertTrue(loaded_model["Main"].HasNodalSolutionStepVariable(KratosMultiphysics.TEMPERATURE))
        self.assertTrue(loaded_model["Other"].HasNodalSolutionStepVariable(KratosMultiphysics.PRESSURE))

        self.assertTrue(loaded_model.HasModelPart("Main.Inlets"))
        self.assertTrue(loaded_model.HasModelPart("Main.Temp"))
        self.assertTrue(1 in loaded_model["Main"].Nodes)
        self.assertTrue(1 in loaded_model["Other"].Nodes)



if __name__ == '__main__':
    KratosUnittest.main()
