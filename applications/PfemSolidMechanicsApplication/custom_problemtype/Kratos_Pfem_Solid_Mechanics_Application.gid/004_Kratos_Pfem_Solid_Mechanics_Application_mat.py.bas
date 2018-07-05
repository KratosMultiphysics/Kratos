
# Importing the Kratos Library
from KratosMultiphysics import **
from KratosMultiphysics.SolidMechanicsApplication import **
from KratosMultiphysics.PfemSolidMechanicsApplication import **
from KratosMultiphysics.ConstitutiveModelsApplication import **
from KratosMultiphysics.UmatApplication import **
#from beam_sections_python_utility import SetProperties

def AssignMaterial(Properties):
*loop materials
# GUI property identifier: Property
*format "%i"
    prop_id = *MatNum;
    prop = Properties[prop_id]
*if(strcmp(MatProp(Type),"FabricModel")==0)
    model = FabricSmallStrainUmatModel()
    mat = SmallStrainAxisymmetric2DLaw(model)
*else
    mat = *MatProp(CONSTITUTIVE_LAW_NAME)()
*endif
    prop.SetValue(CONSTITUTIVE_LAW, mat.Clone())
*end materials
