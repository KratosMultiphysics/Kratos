
from KratosMultiphysics import *
from KratosMultiphysics.SolidMechanicsApplication import *
from KratosMultiphysics.FemToDemApplication import *

def AssignMaterial(Properties):

    prop_id = 1
    prop = Properties[prop_id]
    mat = HyperElasticIsotropicNeoHookeanPlaneStrain2D()
    prop.SetValue(CONSTITUTIVE_LAW, mat.Clone())

    prop_id = 2
    prop = Properties[prop_id]
    mat = HyperElasticIsotropicNeoHookeanPlaneStrain2D()
    prop.SetValue(CONSTITUTIVE_LAW, mat.Clone())

