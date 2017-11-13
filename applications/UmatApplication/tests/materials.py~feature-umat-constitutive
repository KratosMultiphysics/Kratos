
from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.SolidMechanicsApplication import *
from KratosMultiphysics.ConstitutiveLawsApplication import *
#from beam_sections_python_utility import SetProperties
def AssignMaterial(Properties):
    
    prop_id = 1;
    prop = Properties[prop_id]
    mat = LinearElasticPlaneStrain2DLaw()
    prop.SetValue(CONSTITUTIVE_LAW, mat)
    
    prop_id = 2
    prop = Properties[prop_id]
    mat = LinearElastic3DLaw();
    mat = Umat()
    prop.SetValue(CONSTITUTIVE_LAW, mat)
    return;
    
    prop_id = 3;
    prop = Properties[prop_id]
    mat = IsotropicDamageSimoJuPlaneStress2DLaw();
    prop.SetValue(CONSTITUTIVE_LAW, mat.Clone());
    
    prop_id = 4;
    prop = Properties[prop_id]
    mat = LinearElasticPlaneStress2DLaw()
    prop.SetValue(CONSTITUTIVE_LAW, mat.Clone())
    
    prop_id = 5;
    prop = Properties[prop_id]
    mat = HyperElasticPlaneStrain2DLaw()
    prop.SetValue(CONSTITUTIVE_LAW, mat.Clone())
    
    prop_id = 6
    prop = Properties[prop_id]
    mat = HyperElastic3DLaw();
    prop.SetValue(CONSTITUTIVE_LAW, mat.Clone());
 
