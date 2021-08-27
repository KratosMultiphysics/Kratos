import KratosMultiphysics as KM
from KratosMultiphysics.NeuralNetworkApplication.add_layer_process import AddLayerProcess


from importlib import import_module
from tensorflow.keras import layers
import tensorflow.keras as keras

def Factory(settings):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return AddTimeDistributedLayerProcess(settings["parameters"])

class AddTimeDistributedLayerProcess(AddLayerProcess):

    def Add(self, model, hp = None):
        """ Processes to add layers to a network. """
        layer = self.layer_class.Build(hp = hp)
        if isinstance(model, keras.Sequential):
            model.add(layers.TimeDistributed(layer))
            return model
        else:
            return layers.TimeDistributed(layer)(model)

