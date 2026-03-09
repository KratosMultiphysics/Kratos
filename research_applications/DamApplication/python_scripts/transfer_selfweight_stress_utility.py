#importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.StructuralMechanicsApplication import *
# from KratosMultiphysics.ExternalSolversApplication import *
from KratosMultiphysics.DamApplication import *

class TransferSelfweightStressToMainModelPartUtility:

    def __init__(self):

        self.TransferUtility = TransferSelfweightStressUtility()

    def Transfer(self,selfweight_model_part,model_part,domain_size):

        self.TransferUtility.Transfer(selfweight_model_part, model_part, domain_size)

    def TransferInitialStress(self,model_part,domain_size):

        self.TransferUtility.TransferInitialStress(model_part, domain_size)

