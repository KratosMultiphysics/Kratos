// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Anne van de Graaf
//

#include "linear_elastic_law.h"

namespace Kratos
{

bool GeoLinearElasticLaw::RequiresInitializeMaterialResponse()
{
    // Avoid that `ConstitutiveLaw::InitializeMaterialResponseCauchy` throws
    // an exception
    return false;
}

ConstitutiveLaw::StrainMeasure GeoLinearElasticLaw::GetStrainMeasure()
{
    return StrainMeasure_Infinitesimal;
}

ConstitutiveLaw::StressMeasure GeoLinearElasticLaw::GetStressMeasure()
{
    return StressMeasure_Cauchy;
}

void GeoLinearElasticLaw::CalculateMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{
    KRATOS_TRY
    // b.- Get Values to compute the constitutive law:
    const Flags& r_options = rValues.GetOptions();

    Vector& r_strain_vector = rValues.GetStrainVector();

    // NOTE: SINCE THE ELEMENT IS IN SMALL STRAINS WE CAN USE ANY STRAIN MEASURE. HERE EMPLOYING THE CAUCHY_GREEN
    if (r_options.IsNot(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN))
    {
        CalculateCauchyGreenStrain(rValues, r_strain_vector);
    }

    if (r_options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR))
    {
        Matrix& r_constitutive_matrix = rValues.GetConstitutiveMatrix();
        CalculateElasticMatrix(r_constitutive_matrix, rValues);
    }

    if (r_options.Is(ConstitutiveLaw::COMPUTE_STRESS))
    {
        Vector& r_stress_vector = rValues.GetStressVector();
        CalculatePK2Stress(r_strain_vector, r_stress_vector, rValues);
    }

    KRATOS_CATCH("")
}

void GeoLinearElasticLaw::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, ConstitutiveLaw)
}

void GeoLinearElasticLaw::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, ConstitutiveLaw)
}

} // namespace Kratos
