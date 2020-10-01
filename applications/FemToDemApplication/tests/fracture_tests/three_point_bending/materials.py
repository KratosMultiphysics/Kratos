
from __future__ import print_function, absolute_import, division
from KratosMultiphysics import *
from KratosMultiphysics.SolidMechanicsApplication import *
from KratosMultiphysics.FemToDemApplication import *

def AssignMaterial(Properties):

    prop_id = 1
    prop = Properties[prop_id]
    mat = LinearPlaneStress()
    prop.SetValue(CONSTITUTIVE_LAW, mat.Clone())

    prop_id = 2
    prop = Properties[prop_id]
    mat = LinearPlaneStress()
    prop.SetValue(CONSTITUTIVE_LAW, mat.Clone())

