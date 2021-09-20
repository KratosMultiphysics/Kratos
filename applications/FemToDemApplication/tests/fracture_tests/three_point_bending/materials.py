
from KratosMultiphysics import *
from KratosMultiphysics.FemToDemApplication import *

def AssignMaterial(Properties):

    prop_id = 1
    prop = Properties[prop_id]
    mat = LinearPlaneStressFEMDEM()
    prop.SetValue(CONSTITUTIVE_LAW, mat.Clone())

    prop_id = 2
    prop = Properties[prop_id]
    mat = LinearPlaneStressFEMDEM()
    prop.SetValue(CONSTITUTIVE_LAW, mat.Clone())

