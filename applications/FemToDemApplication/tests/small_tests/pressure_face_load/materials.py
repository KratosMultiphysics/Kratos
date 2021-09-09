
from __future__ import print_function, absolute_import, division
from KratosMultiphysics import *
from KratosMultiphysics.FemToDemApplication import *

def AssignMaterial(Properties):

    prop_id = 1
    prop = Properties[prop_id]
    mat = ElasticIsotropic3DFEMDEM()
    prop.SetValue(CONSTITUTIVE_LAW, mat.Clone())

