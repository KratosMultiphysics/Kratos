# importing the Kratos Library
import KratosMultiphysics
import materials

class ConstitutiveLawUtility:
    #

    def __init__(self, model_part, dimension):

        self.model_part = model_part
        self.dimension  = dimension

    #
    def Initialize(self):
        self.SetConstitutiveLaw()

    #
    def SetConstitutiveLaw(self):
        materials.AssignMaterial(self.model_part.Properties)



    #
