// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Alejandro Cornejo
//
// System includes

// External includes

// Project includes
#include "timoshenko_beam_elastic_constitutive_law.h"
#include "structural_mechanics_application_variables.h"

namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

TimoshenkoBeamElasticConstitutiveLaw::TimoshenkoBeamElasticConstitutiveLaw()
    : BeamConstitutiveLaw()
{
}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

TimoshenkoBeamElasticConstitutiveLaw::TimoshenkoBeamElasticConstitutiveLaw(const TimoshenkoBeamElasticConstitutiveLaw& rOther)
    : BeamConstitutiveLaw(rOther)
{
}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer TimoshenkoBeamElasticConstitutiveLaw::Clone() const
{
    TimoshenkoBeamElasticConstitutiveLaw::Pointer p_clone(new TimoshenkoBeamElasticConstitutiveLaw(*this));
    return p_clone;
}

//*************************CONSTITUTIVE LAW GENERAL FEATURES *************************
//************************************************************************************

void TimoshenkoBeamElasticConstitutiveLaw::GetLawFeatures(Features& rFeatures)
{
    rFeatures.mStrainSize = 3;
    rFeatures.mSpaceDimension = 2;
}

//************************************************************************************
//************************************************************************************

int TimoshenkoBeamElasticConstitutiveLaw::Check(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const ProcessInfo& rCurrentProcessInfo
) const
{
    KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(YOUNG_MODULUS)) << "YOUNG_MODULUS is not defined in the properties" << std::endl;
    KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(POISSON_RATIO)) << "POISSON_RATIO is not defined in the properties" << std::endl;
    KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(CROSS_AREA))    << "CROSS_AREA is not defined in the properties"    << std::endl;
    KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(DENSITY))       << "DENSITY is not defined in the properties"       << std::endl;
    KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(I33))           << "I33 is not defined in the properties"           << std::endl;
    return 0;
}

//************************************************************************************
//************************************************************************************

void TimoshenkoBeamElasticConstitutiveLaw::CalculateMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{
    CalculateMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

void TimoshenkoBeamElasticConstitutiveLaw::CalculateMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues)
{
    CalculateMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

void TimoshenkoBeamElasticConstitutiveLaw::CalculateMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues)
{
    CalculateMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

void TimoshenkoBeamElasticConstitutiveLaw::CalculateMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
    auto &r_material_properties = rValues.GetMaterialProperties();
    auto &r_strain_vector = rValues.GetStrainVector();

    const double axial_strain = r_strain_vector[0]; // E_l
    const double curvature    = r_strain_vector[1]; // Kappa
    const double shear_strain = r_strain_vector[2]; // Gamma_xy

    const double E  = r_material_properties[YOUNG_MODULUS];
    const double nu = r_material_properties[POISSON_RATIO];
    const double A  = r_material_properties[CROSS_AREA];
    const double I  = r_material_properties[IZ];
}

} // Namespace Kratos
