from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
import KratosMultiphysics
import materials
KratosMultiphysics.CheckForPreviousImport()


class ConstitutiveLawUtility:
    #

    def __init__(self, model_part, domain_size):

        self.model_part  = model_part
        self.domain_size = domain_size

    #
    def Initialize(self):
        self.SetConstitutiveLaw()

    #
    def SetConstitutiveLaw(self):
        materials.AssignMaterial(self.model_part.Properties)



    #
