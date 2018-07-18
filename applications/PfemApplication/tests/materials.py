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
    #multiple solid material definitions:
    #V-P element
    #mat = KratosMaterial.StrainRatePlaneStrain2DLaw(KratosMaterial.IncompressibleHypoElasticModel())
    #V element
    #mat = KratosMaterial.StrainRatePlaneStrain2DLaw(KratosMaterial.IsochoricHypoElasticModel())
    #mat = KratosMaterial.StrainRatePlaneStrain2DLaw(KratosMaterial.HypoElasticModel())
    mat = KratosMaterial.LargeStrainPlaneStrain2DLaw(KratosMaterial.SaintVenantKirchhoffModel())
    prop.SetValue(KratosMultiphysics.CONSTITUTIVE_LAW, mat.Clone())

    prop_id = 3
    prop = Properties[prop_id]
    mat = KratosMaterial.Newtonian3DLaw()
    prop.SetValue(KratosMultiphysics.CONSTITUTIVE_LAW, mat.Clone())
