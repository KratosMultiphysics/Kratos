import KratosMultiphysics as KM
import KratosMultiphysics.NeuralNetworkApplication as NeuralNetworkApplication
import KratosMultiphysics.KratosUnittest as KratosUnittest
import tensorflow as tf
from KratosMultiphysics.NeuralNetworkApplication.initiate_network_process import InitiateNetworkProcess
from KratosMultiphysics.NeuralNetworkApplication.add_layer_process import AddLayerProcess

class TestInitialLayers(KratosUnittest.TestCase):
    def _execute_initialization_layer_test(self,layer):
        parameters = KM.Parameters()
        parameters.AddEmptyValue("layer_type")
        parameters["layer_type"].SetString(layer + "_layer")
        parameters.AddEmptyValue("layer_parameters")
        parameters["layer_parameters"].AddEmptyValue("data_input")
        parameters["layer_parameters"]["data_input"].SetString("data/training_in_raw.dat")
        layer_process = InitiateNetworkProcess(parameters)
        layer = layer_process.Initialize()
        self.assertEqual(layer.shape.as_list(), [None, 1])

    def test_Input(self):
        with KratosUnittest.WorkFolderScope(".", __file__):
            self._execute_initialization_layer_test(layer = 'input')

class TestInnerLayers(KratosUnittest.TestCase):
    def _execute_inner_layer_test(self,layer,neurons):
        model = tf.keras.Sequential()
        model.add(tf.keras.Input(shape=(1,)))
        parameters = KM.Parameters()
        parameters.AddEmptyValue("layer_type")
        parameters["layer_type"].SetString(layer + "_layer")
        parameters.AddEmptyValue("layer_parameters")
        parameters["layer_parameters"].AddEmptyValue("units")
        parameters["layer_parameters"]["units"].SetInt(neurons)
        layer_process = AddLayerProcess(parameters)
        model = layer_process.Add(model)
        self.assertEqual(model.output_shape, (None, neurons))

    def test_DenseLayer(self):
        self._execute_inner_layer_test(layer = 'dense', neurons = 4)

class TestBlockLayers(KratosUnittest.TestCase):
    def _execute_block_layer_test(self,layer,architecture):
        model = tf.keras.Sequential()
        model.add(tf.keras.Input(shape=(1,)))
        parameters = KM.Parameters()
        parameters.AddEmptyValue("layer_type")
        parameters["layer_type"].SetString(layer + "_layer")
        parameters.AddEmptyValue("layer_parameters")
        parameters["layer_parameters"].AddEmptyValue("number_layers")
        parameters["layer_parameters"]["number_layers"].SetInt(len(architecture))
        parameters["layer_parameters"].AddEmptyValue("neurons_per_layer")
        parameters["layer_parameters"]["neurons_per_layer"].SetVector(architecture)
        layer_process = AddLayerProcess(parameters)
        model = layer_process.Add(model)
        self.assertEqual(model.output_shape, (None, architecture[-1]))
    def test_BlockDenseLayer(self):
        self._execute_block_layer_test(layer = 'block_dense', architecture = [16,4,2])

if __name__ == '__main__':
    KratosUnittest.main()