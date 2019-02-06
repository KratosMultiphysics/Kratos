
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
*if(strcmp(MatProp(CONSTITUTIVE_LAW_NAME),"FabricModelAxisym2DLaw")==0)
    mat = SmallStrainAxisymmetric2DLaw(model)
*else
    mat = SmallStrainPlaneStrain2DLaw(model)
*endif
*elseif(strcmp(MatProp(Type),"GensNovaPlasticity")==0)
*if(strcmp(MatProp(CONSTITUTIVE_LAW_NAME),"GensNovaPlasticAxisym2DLaw")==0)
    model = GensNovaModel()
    mat = LargeStrainAxisymmetric2DLaw(model)
*elseif(strcmp(MatProp(CONSTITUTIVE_LAW_NAME),"V2GensNovaPlasticAxisym2DLaw")==0)
    model = V2GensNovaModel()
    mat = LargeStrainAxisymmetric2DLaw(model)
*elseif(strcmp(MatProp(CONSTITUTIVE_LAW_NAME),"NonlocalV2GensNovaPlasticAxisym2DLaw")==0)
    model = NonlocalV2GensNovaModel()
    mat = LargeStrainAxisymmetric2DLaw(model)
*elseif(strcmp(MatProp(CONSTITUTIVE_LAW_NAME),"V2GensNovaPlasticPlaneStrain2DLaw")==0)
    model = V2GensNovaModel()
    mat = LargeStrainPlaneStrain2DLaw(model)
*elseif(strcmp(MatProp(CONSTITUTIVE_LAW_NAME),"NonlocalV2GensNovaPlasticPlaneStrain2DLaw")==0)
    model = NonlocalV2GensNovaModel()
    mat = LargeStrainPlaneStrain2DLaw(model)
*else
    model = GensNovaModel()
    mat = LargeStrainPlaneStrain2DLaw(model)
*endif
*else
    mat = *MatProp(CONSTITUTIVE_LAW_NAME)()
*endif
    prop.SetValue(CONSTITUTIVE_LAW, mat.Clone())
*end materials
