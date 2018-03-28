
# Importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.SolidMechanicsApplication import *
from KratosMultiphysics.PfemSolidMechanicsApplication import *
#from beam_sections_python_utility import SetProperties

def AssignMaterial(Properties):
# GUI property identifier: Property
    prop_id = 1;
    prop = Properties[prop_id]
    mat = HenckyTrescaPlasticPlaneStrain2DLaw()
    mat = HyperElasticPlaneStrain2DLaw()
    mat = LinearElasticPlaneStrain2DLaw()
    prop.SetValue(CONSTITUTIVE_LAW, mat.Clone())
