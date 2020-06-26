from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.StructuralMechanicsApplication import *
# from KratosMultiphysics.ExternalSolversApplication import *
from KratosMultiphysics.DamApplication import *

class GlobalJoinStressUtility:

    def __init__(self,model_part):

        # This is an example of input paramerter.
        params = Parameters("{}")
        params.AddEmptyValue("pmid0_0").SetDouble(247.67)
        params.AddEmptyValue("pmid0_1").SetDouble(-368.92)
        params.AddEmptyValue("pmid0_2").SetDouble(242.6)
        params.AddEmptyValue("pmid1_0").SetDouble(247.15)
        params.AddEmptyValue("pmid1_1").SetDouble(-385.15)
        params.AddEmptyValue("pmid1_2").SetDouble(242.6)
        params.AddEmptyValue("pmid2_0").SetDouble(313.88)
        params.AddEmptyValue("pmid2_1").SetDouble(-373.26)
        params.AddEmptyValue("pmid2_2").SetDouble(242.6)

        self.JointStressUtility = GlobalJointStressUtility(model_part, params)

    def ComputingGlobalStress(self):

        self.JointStressUtility.ComputingGlobalStress()
