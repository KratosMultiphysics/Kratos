
from __future__ import print_function, absolute_import, division
from KratosMultiphysics import *
from KratosMultiphysics.SolidMechanicsApplication import *
from KratosMultiphysics.FemToDemApplication import *

def AssignMaterial(Properties):

    prop_id = 1
    prop = Properties[prop_id]
    mat = HyperElasticIsotropicNeoHookeanPlaneStrain2D()
    prop.SetValue(CONSTITUTIVE_LAW, mat.Clone())

