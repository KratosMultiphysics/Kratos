
# Importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.SolidMechanicsApplication import *
from KratosMultiphysics.LagrangianMPMApplication import *
#from beam_sections_python_utility import SetProperties
def AssignMaterial(Properties):
# GUI property identifier: Property1
# GUI material identifier: Steel_AISI1059
    prop_id = 1;
    prop = Properties[prop_id]
    #mat = HyperElasticPlasticJ2PlaneStrain2DLaw();#HyperElasticPlasticJ2Axisym2DLaw();#
    mat = LinearElasticPlaneStress2DLaw()
    prop.SetValue(CONSTITUTIVE_LAW, mat.Clone());
    print(prop)
