
class ElementSelectionStrategy():
    def __init__(self):
        """ The required tolerances are entered here. Default values to be provided
        """

    def SetUp(self):
        """ This method treats the data passed to the object, in order to obtain a suitable format
        """

    def Initialize(self):
        """ This method performs the initialization of the variables to be used
        """

    def Calculate(self):
        """ This method obtains a set of elements and conditions to be used for hyper-reduction
        """

    def WriteSelectedElements(self):
        """ This method writes to the RomParameters.json file the selected elements and conditions to be used in the hyper-reduction
        simulation
        """

    def Run(self):
        self.Initialize()
        self.Calculate()
        self.WriteSelectedElements()



