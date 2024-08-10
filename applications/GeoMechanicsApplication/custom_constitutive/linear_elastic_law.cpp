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
#include "includes/mat_variables.h"

namespace Kratos
{

bool GeoLinearElasticLaw::RequiresInitializeMaterialResponse()
{
    // Avoid that `ConstitutiveLaw::InitializeMaterialResponse` throws
    // an exception
    return false;
}

bool GeoLinearElasticLaw::RequiresFinalizeMaterialResponse()
{
    // Avoid that `ConstitutiveLaw::FinalizeMaterialResponse` throws
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

// NOTE: Since we are in the hypothesis of small strains we can use the same function for everything
void GeoLinearElasticLaw::CalculateMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues)
{
    CalculateMaterialResponsePK2(rValues);
}

void GeoLinearElasticLaw::CalculateMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{
    KRATOS_TRY

    const Flags& r_options = rValues.GetOptions();

    Vector& r_strain_vector = rValues.GetStrainVector();

    // NOTE: SINCE THE ELEMENT IS IN SMALL STRAINS WE CAN USE ANY STRAIN MEASURE. HERE EMPLOYING THE CAUCHY_GREEN
    if (r_options.IsNot(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN)) {
        CalculateCauchyGreenStrain(rValues, r_strain_vector);
    }

    if (r_options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
        Matrix& r_constitutive_matrix = rValues.GetConstitutiveMatrix();
        CalculateElasticMatrix(r_constitutive_matrix, rValues);
    }

    if (r_options.Is(ConstitutiveLaw::COMPUTE_STRESS)) {
        Vector& r_stress_vector = rValues.GetStressVector();
        CalculatePK2Stress(r_strain_vector, r_stress_vector, rValues);
    }

    KRATOS_CATCH("")
}

void GeoLinearElasticLaw::CalculateMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues)
{
    CalculateMaterialResponsePK2(rValues);
}

void GeoLinearElasticLaw::CalculateMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
    CalculateMaterialResponsePK2(rValues);
}

double& GeoLinearElasticLaw::CalculateValue(ConstitutiveLaw::Parameters& rParameterValues,
                                            const Variable<double>&      rThisVariable,
                                            double&                      rValue)
{
    if (rThisVariable == STRAIN_ENERGY) {
        Vector& r_strain_vector = rParameterValues.GetStrainVector();
        Vector& r_stress_vector = rParameterValues.GetStressVector();
        this->CalculateCauchyGreenStrain(rParameterValues, r_strain_vector);
        this->CalculatePK2Stress(r_strain_vector, r_stress_vector, rParameterValues);

        rValue = 0.5 * inner_prod(r_strain_vector, r_stress_vector); // Strain energy = 0.5*E:C:E
    }

    return rValue;
}

Vector& GeoLinearElasticLaw::CalculateValue(ConstitutiveLaw::Parameters& rParameterValues,
                                            const Variable<Vector>&      rThisVariable,
                                            Vector&                      rValue)
{
    if (rThisVariable == STRAIN || rThisVariable == GREEN_LAGRANGE_STRAIN_VECTOR ||
        rThisVariable == ALMANSI_STRAIN_VECTOR) {
        this->CalculateCauchyGreenStrain(rParameterValues, rValue);
    } else if (rThisVariable == STRESSES || rThisVariable == CAUCHY_STRESS_VECTOR ||
               rThisVariable == KIRCHHOFF_STRESS_VECTOR || rThisVariable == PK2_STRESS_VECTOR) {
        // Get Values to compute the constitutive law:
        Flags& r_flags = rParameterValues.GetOptions();

        // Previous flags saved
        const bool flag_strain       = r_flags.Is(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);
        const bool flag_const_tensor = r_flags.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
        const bool flag_stress       = r_flags.Is(ConstitutiveLaw::COMPUTE_STRESS);

        r_flags.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, false);
        r_flags.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);
        r_flags.Set(ConstitutiveLaw::COMPUTE_STRESS, true);

        // We compute the stress
        CalculateMaterialResponseCauchy(rParameterValues);
        rValue = rParameterValues.GetStressVector();

        // Previous flags restored
        r_flags.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, flag_strain);
        r_flags.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, flag_const_tensor);
        r_flags.Set(ConstitutiveLaw::COMPUTE_STRESS, flag_stress);
    }

    return rValue;
}

Matrix& GeoLinearElasticLaw::CalculateValue(ConstitutiveLaw::Parameters& rParameterValues,
                                            const Variable<Matrix>&      rThisVariable,
                                            Matrix&                      rValue)
{
    if (rThisVariable == CONSTITUTIVE_MATRIX || rThisVariable == CONSTITUTIVE_MATRIX_PK2 ||
        rThisVariable == CONSTITUTIVE_MATRIX_KIRCHHOFF) {
        CalculateElasticMatrix(rValue, rParameterValues);
    }

    return rValue;
}

void GeoLinearElasticLaw::SetValue(const Variable<double>&, const double&, const ProcessInfo&)
{
    // Since `ConstitutiveLaw::SetValue` always throws an exception, an override
    // is required when that is undesired behavior
}

void GeoLinearElasticLaw::SetValue(const Variable<Vector>&, const Vector&, const ProcessInfo&)
{
    // Since `ConstitutiveLaw::SetValue` always throws an exception, an override
    // is required when that is undesired behavior
}

int GeoLinearElasticLaw::Check(const Properties&   rMaterialProperties,
                               const GeometryType& rElementGeometry,
                               const ProcessInfo&  rCurrentProcessInfo) const
{
    KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(YOUNG_MODULUS))
        << "YOUNG_MODULUS is not available in material parameters" << std::endl;
    KRATOS_ERROR_IF(rMaterialProperties[YOUNG_MODULUS] <= 0.0)
        << "YOUNG_MODULUS is invalid value " << std::endl;

    KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(POISSON_RATIO))
        << "POISSON_RATIO is not available in material parameters" << std::endl;

    const double& nu = rMaterialProperties[POISSON_RATIO];
    KRATOS_ERROR_IF((nu > 0.499 && nu < 0.501) || (nu < -0.999 && nu > -1.01))
        << "POISSON_RATIO has invalid value " << std::endl;

    return 0;
}

void GeoLinearElasticLaw::SetConsiderDiagonalEntriesOnlyAndNoShear(bool Whether)
{
    mConsiderDiagonalEntriesOnlyAndNoShear = Whether;
}

bool GeoLinearElasticLaw::GetConsiderDiagonalEntriesOnlyAndNoShear() const
{
    return mConsiderDiagonalEntriesOnlyAndNoShear;
}

void GeoLinearElasticLaw::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, ConstitutiveLaw)
    // Since `mConsiderDiagonalEntriesOnlyAndNoShear` is supposed to be changed only during the
    // execution of a stage that includes the K0 procedure, there is no need to serialize it
}

void GeoLinearElasticLaw::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, ConstitutiveLaw)
    // Since `mConsiderDiagonalEntriesOnlyAndNoShear` is supposed to be changed only during the
    // execution of a stage that includes the K0 procedure, there is no need to serialize it
}

} // namespace Kratos
