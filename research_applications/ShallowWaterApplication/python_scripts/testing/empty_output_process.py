# Importing the KratosMultiphysics library
import KratosMultiphysics as KM

def Factory(settings, Model):
    return EmptyOutputProcess(Model, settings["Parameters"])

class EmptyOutputProcess(KM.OutputProcess):
    def __init__(self, model, parameters):
        ''' EmptyOutputProcess constructor.

        This is a temporal file and it must be deleted once there is an output_process base class in the core.
        The main reason of existing is the need to test VisualizationMeshProcess with JsonOutputProcess and FromJsonCheckResultProcess.
        '''
        KM.OutputProcess.__init__(self)

    def PrintOutput(self):
        pass
