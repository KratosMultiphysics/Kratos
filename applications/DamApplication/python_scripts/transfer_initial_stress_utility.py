#importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.StructuralMechanicsApplication import *
# from KratosMultiphysics.ExternalSolversApplication import *
from KratosMultiphysics.DamApplication import *

class TransferInitialStressToMainModelPartUtility:

    def __init__(self):

        self.TransferUtility = TransferInitialStressUtility()

    def TransferInitialStress(self,model_part,domain_size):

        self.TransferUtility.TransferInitialStress(model_part, domain_size)
