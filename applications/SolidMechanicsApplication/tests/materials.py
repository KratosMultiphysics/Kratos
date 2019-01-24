# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.ConstitutiveModelsApplication as KratosMaterial

def AssignMaterial(Properties):

    prop_id = 1
    prop = Properties[prop_id]
    mat = KratosMaterial.SmallStrainPlaneStrain2DLaw(KratosMaterial.LinearElasticModel())
    prop.SetValue(KratosMultiphysics.CONSTITUTIVE_LAW, mat.Clone())

    prop_id = 2
    prop = Properties[prop_id]
    mat = KratosMaterial.SmallStrain3DLaw(KratosMaterial.LinearElasticModel())
    prop.SetValue(KratosMultiphysics.CONSTITUTIVE_LAW, mat.Clone())

    prop_id = 3
    prop = Properties[prop_id]
    mat = KratosMaterial.SmallStrainPlaneStress2DLaw(KratosMaterial.SimoJuExponentialDamageModel())
    prop.SetValue(KratosMultiphysics.CONSTITUTIVE_LAW, mat.Clone())

    prop_id = 4
    prop = Properties[prop_id]
    mat = KratosMaterial.SmallStrainPlaneStress2DLaw(KratosMaterial.LinearElasticModel())
    prop.SetValue(KratosMultiphysics.CONSTITUTIVE_LAW, mat.Clone())

    prop_id = 5
    prop = Properties[prop_id]
    mat = KratosMaterial.LargeStrainPlaneStrain2DLaw(KratosMaterial.SaintVenantKirchhoffModel())
    prop.SetValue(KratosMultiphysics.CONSTITUTIVE_LAW, mat.Clone())

    prop_id = 6
    prop = Properties[prop_id]
    mat = KratosMaterial.LargeStrain3DLaw(KratosMaterial.SaintVenantKirchhoffModel())
    prop.SetValue(KratosMultiphysics.CONSTITUTIVE_LAW, mat.Clone())

