
from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# Importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.ConstitutiveModelsApplication import *
#from beam_sections_python_utility import SetProperties
#from beam_sections_python_utility import SetMaterialProperties

def AssignMaterial(Properties):
    # material for solid material

    prop_id = 1;
    prop = Properties[prop_id]
    young_modulus = prop.GetValue(YOUNG_MODULUS)
    poisson_ratio = prop.GetValue(POISSON_RATIO)
    c10 = young_modulus*0.25/(1.0+poisson_ratio)
    print(c10)
    prop.SetValue( C10, c10)
    mat = LargeStrainAxisymmetric2DLaw(NeoHookeanModel())
    prop.SetValue(CONSTITUTIVE_LAW, mat.Clone())
        
        
