from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.StructuralMechanicsApplication import *
from KratosMultiphysics.PoromechanicsApplication import *
from KratosMultiphysics.DamApplication import *

def SetConstitutiveLaw(model_part):

    for prop in model_part.Properties:
        if prop.Has(CONSTITUTIVE_LAW_NAME):
            ConstitutiveLawName=prop.GetValue(CONSTITUTIVE_LAW_NAME)
            if (ConstitutiveLawName=="Newtonian"):
                print("Newtonian Law is not used, we just need its parameters")
            else:
                mat = globals().get(ConstitutiveLawName)()
                prop.SetValue(CONSTITUTIVE_LAW, mat.Clone())
        else:
            print("Property",prop.Id,"has no CONSTITUTIVE_LAW_NAME")
