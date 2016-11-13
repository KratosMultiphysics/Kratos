from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.SolidMechanicsApplication import *
from KratosMultiphysics.DamApplication import *
from KratosMultiphysics.PoromechanicsApplication import *
CheckForPreviousImport()


def SetConstitutiveLaw(model_part):
    for prop in model_part.Properties:
        ConstitutiveLawName=prop.GetValue(CONSTITUTIVE_LAW_NAME)
        if(ConstitutiveLawName == "LinearElastic2DPlaneStress"):
            prop.SetValue(CONSTITUTIVE_LAW, LinearElasticPlaneStress2DLaw())
        elif(ConstitutiveLawName == "LinearElastic2DPlaneStrain"):
            prop.SetValue(CONSTITUTIVE_LAW, LinearElasticPlaneStrain2DLaw())
        elif(ConstitutiveLawName == "LinearElastic3D"):
            prop.SetValue(CONSTITUTIVE_LAW, LinearElastic3DLaw())
        elif(ConstitutiveLawName == "ThermalLinearElastic2DPlaneStress"):
            prop.SetValue(CONSTITUTIVE_LAW, ThermalLinearElastic2DPlaneStress())
        elif(ConstitutiveLawName == "ThermalLinearElastic2DPlaneStrain"):
            prop.SetValue(CONSTITUTIVE_LAW, ThermalLinearElastic2DPlaneStrain())
        elif(ConstitutiveLawName == "ThermalLinearElastic3D"):
            prop.SetValue(CONSTITUTIVE_LAW, ThermalLinearElastic3DLaw())
        elif(ConstitutiveLawName == "BilinearCohesive2DLaw"):
            prop.SetValue(CONSTITUTIVE_LAW, BilinearCohesive2DLaw())
        elif(ConstitutiveLawName == "BilinearCohesive3DLaw"):
            prop.SetValue(CONSTITUTIVE_LAW, BilinearCohesive3DLaw())
        elif(ConstitutiveLawName == "SimoJuLocalDamage3DLaw"):
            prop.SetValue(CONSTITUTIVE_LAW, SimoJuLocalDamage3DLaw())
        elif(ConstitutiveLawName == "SimoJuLocalDamagePlaneStrain2DLaw"):
            prop.SetValue(CONSTITUTIVE_LAW, SimoJuLocalDamagePlaneStrain2DLaw())
        elif(ConstitutiveLawName == "SimoJuLocalDamagePlaneStress2DLaw"):
            prop.SetValue(CONSTITUTIVE_LAW, SimoJuLocalDamagePlaneStress2DLaw())

