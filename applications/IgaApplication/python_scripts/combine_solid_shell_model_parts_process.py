# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.IgaApplication as IGA


def Factory(settings, model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    return CombineSolidShellModelPartsProcess(model, settings["Parameters"])

class CombineSolidShellModelPartsProcess(KratosMultiphysics.Process):
    """Only the member variables listed below should be accessed directly.

    Public member variables:
    model -- the container of the different model parts.
    params -- Kratos parameters containing the settings.
    """
    def __init__(self,model, params):
        """ The default constructor of the class

        Keyword arguments:
        self -- It signifies an instance of a class.
        model -- the container of the different model parts.
        params -- Kratos parameters containing solver settings.
        """
        #print("Minas: Constructor of CombineSolidShellModelPartsProcess!!!")

        KratosMultiphysics.Process.__init__(self)

        ## Overwrite the default settings with user-provided parameters
        self.CombinedModel = model
        self.CombineSolidShellModelPartsProcess =  IGA.CombineSolidShellModelPartsProcess(self.CombinedModel)

    def ExecuteInitialize(self):
        # Get the model parts which divide the problem
        #print("Minas: Executing CombineSolidShellModelPartsProcess")
        self.CombineSolidShellModelPartsProcess.ExecuteInitialize()