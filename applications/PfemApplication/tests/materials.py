# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.ConstitutiveModelsApplication as KratosMaterial

def AssignMaterial(Properties):
    
    prop_id = 1
    prop = Properties[prop_id]
    mat = KratosMaterial.NewtonianPlaneStrain2DLaw()
    prop.SetValue(KratosMultiphysics.CONSTITUTIVE_LAW, mat.Clone())
    
    prop_id = 2
    prop = Properties[prop_id]
    mat = KratosMaterial.Newtonian3DLaw()
    prop.SetValue(KratosMultiphysics.CONSTITUTIVE_LAW, mat.Clone())
    
 
