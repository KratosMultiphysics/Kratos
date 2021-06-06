# Importing the KratosMultiphysics library
import KratosMultiphysics as KM

def Factory(settings, Model):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return EmptyOutputProcess(Model, settings["Parameters"])

class EmptyOutputProcess(KM.Process):
    def __init__(self, model, parameters):
        ''' EmptyOutputProcess constructor.

        This is a temporal file and it must be deleted once there is an output_process base class in the core.
        The main reason of existing is the need to test VisualizationMeshProcess with JsonOutputProcess and FromJsonCheckResultProcess.
        '''
        KM.Process.__init__(self)

    def IsOutputStep(self):
        return True

    def PrintOutput(self):
        pass
