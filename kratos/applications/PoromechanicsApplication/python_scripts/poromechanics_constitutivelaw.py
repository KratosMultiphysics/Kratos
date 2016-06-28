from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.SolidMechanicsApplication import *
from KratosMultiphysics.PoromechanicsApplication import *

def SetConstitutiveLaw(model_part):
    
    for prop in model_part.Properties:
        ConstitutiveLawName=prop.GetValue(CONSTITUTIVE_LAW_NAME)
        mat = globals().get(ConstitutiveLawName)()
        prop.SetValue(CONSTITUTIVE_LAW, mat.Clone())