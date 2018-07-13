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
    #mat = KratosMaterial.StrainRatePlaneStrain2DLaw(KratosMaterial.IncompressibleHypoElasticModel())
    mat = KratosMaterial.LargeStrainPlaneStrain2DLaw(KratosMaterial.SaintVenantKirchhoffModel())    
    prop.SetValue(KratosMultiphysics.CONSTITUTIVE_LAW, mat.Clone())
