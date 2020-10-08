
from KratosMultiphysics import *
from KratosMultiphysics.SolidMechanicsApplication import *
from KratosMultiphysics.FemToDemApplication import *

def AssignMaterial(Properties):

    prop_id = 1
    prop = Properties[prop_id]
    mat = ElasticIsotropic3D()
    prop.SetValue(CONSTITUTIVE_LAW, mat.Clone())

