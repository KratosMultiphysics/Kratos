#importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.SolidMechanicsApplication import *
CheckForPreviousImport()

class ConstitutiveLawUtility:
    #######################################################################
    def __init__(self,model_part,domain_size):

        self.model_part  = model_part
        self.domain_size = domain_size 
        

    #######################################################################
    def Initialize(self):
        self.SetConstitutiveLaw();

    #######################################################################   
    def SetConstitutiveLaw(self):
        from materials import *
        AssignMaterial(self.model_part.Properties)



    #######################################################################   


