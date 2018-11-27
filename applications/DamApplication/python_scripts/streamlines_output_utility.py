from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.SolidMechanicsApplication import *
# from KratosMultiphysics.ExternalSolversApplication import *
from KratosMultiphysics.DamApplication import *

class StreamlinesOutputUtility:

    def __init__(self,domain_size):

        # Construct the utility
        self.domain_size = domain_size
        if domain_size==3:
            self.OutputUtility = StreamlinesOutput3DUtilities()

        else:
            print("**************** Must be implemented ****************")


    def ComputeOutputStep(self, main_model_part ,domain_size):

        self.OutputUtility.ComputeOutputStep( main_model_part ,domain_size)
