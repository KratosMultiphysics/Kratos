
# Importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.StructuralMechanicsApplication import *
#from beam_sections_python_utility import SetProperties
def AssignMaterial(Properties):
    
    prop_id = 1;
    prop = Properties[prop_id]
    mat = LinearElasticPlaneStrain2DLaw()
    prop.SetValue(CONSTITUTIVE_LAW, mat.Clone())
    
    prop_id = 2
    prop = Properties[prop_id]
    mat = LinearElasticIsotropic3DLaw();
    prop.SetValue(CONSTITUTIVE_LAW, mat.Clone());
    
    prop_id = 3;
    prop = Properties[prop_id]
    mat = LinearElasticPlaneStrain2DLaw() ##IsotropicDamageSimoJuPlaneStress2DLaw();
    prop.SetValue(CONSTITUTIVE_LAW, mat.Clone());
    
    prop_id = 4;
    prop = Properties[prop_id]
    mat = LinearElasticPlaneStress2DLaw()
    prop.SetValue(CONSTITUTIVE_LAW, mat.Clone())
    

