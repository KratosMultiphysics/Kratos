
# Importing the Kratos Library
from KratosMultiphysics import *
import KratosMultiphysics.ExternalSolversApplication
from KratosMultiphysics.SolidMechanicsApplication import *
from  KratosMultiphysics.ConstitutiveModelsApplication import *
from KratosMultiphysics.UmatApplication import *
import KratosMultiphysics.PfemApplication
import KratosMultiphysics.ContactMechanicsApplication
from KratosMultiphysics.PfemSolidMechanicsApplication import *
#from beam_sections_python_utility import SetProperties

def AssignMaterial(Properties):
# GUI property identifier: Property
    prop_id = 1;
    prop = Properties[prop_id]
    model = FabricSmallStrainUmatModel()
    mat = SmallStrainAxisymmetric2DLaw(model)
    prop.SetValue(CONSTITUTIVE_LAW, mat.Clone())
