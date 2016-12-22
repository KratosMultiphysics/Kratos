
from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# Importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.SolidMechanicsApplication import *
#from beam_sections_python_utility import SetProperties
#from beam_sections_python_utility import SetMaterialProperties

def AssignMaterial(Properties):
    # material for solid material

    prop_id = 1;
    prop = Properties[prop_id]
    mat = LinearElastic3DLaw()
    prop.SetValue(CONSTITUTIVE_LAW, mat.Clone())
        
