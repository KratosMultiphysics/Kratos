import KratosMultiphysics as KM

def Factory(settings):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a parameters object, encapsulating a json string")
    return NeuralNetworkProcess(settings["parameters"])


class NeuralNetworkProcess(KM.Process):
    """
    This class is the parent for the classes that deal with the training of the neural network.

    It is similar to the normal Kratos Process but it does not receive a model, only settings.
    """
    def __init__(self):
        super().__init__()

    def ExecuteInitialize(self):
        """ Processes to act on the initialization. """
        pass

    def ExecuteFinalize(self):
        """ Processes to act on the finalization. """
        pass

    def ExecuteBeforeTraining(self):
        """ Processes to act just before the training. """
        pass
    
    def ExecuteTraining(self, loss_function, optimizer, callbacks_list, metrics_list):
        """ Processes to act directly during the training step. """
        pass

    def ExecuteTuning(self, hypermodel):
        """ Process to act directly during the tuning step. """
        pass
    
    def Preprocess(self, data_in, data_out):
        """ Process to preprocess the data. """
        pass

    def Invert(self, data_in, data_out):
        """ Process to invert the transformations of the data. """
        pass

    def TransformPredictions(self, processes):
        """ Process to transform the predictions."""
        pass

    def Initialize(self):
        """ Process to initialize a network. """
        pass

    def Add(self,model, hp = None):
        """ Process for adding layers to the network. """
        return model

    def Save(self,model):
        """ Process for saving a network. """
        pass

    def Compile(self, loss_function, optimizer, hp = None):
        """ Process for compiling a network. """
        pass

    def Callback(self):
        """ Process for setting up a callback. """
        pass

    def CompileMetric(self):
        """ Process for compiling a metric. """
        pass

    def Plot(self):
        """ Process for plotting."""
        pass