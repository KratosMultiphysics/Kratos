# Importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.SolidMechanicsApplication import *
def AssignMaterial(Properties):
    prop_id = 1;
    prop = Properties[1]
    mat = LinearElasticPlaneStrain2DLaw();
    prop.SetValue(CONSTITUTIVE_LAW, mat.Clone());
