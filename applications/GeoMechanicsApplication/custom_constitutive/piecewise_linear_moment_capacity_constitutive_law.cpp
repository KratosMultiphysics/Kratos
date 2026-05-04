// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//
//  Main authors:    Gennady Markelov

// Project includes
#include "piecewise_linear_moment_capacity_constitutive_law.h"
#include "custom_utilities/check_utilities.hpp"
#include "geo_mechanics_application_variables.h"
#include "includes/properties.h"

#include <algorithm>

namespace Kratos
{
using namespace std::string_literals;

ConstitutiveLaw::Pointer PiecewiseLinearMomentCapacityConstitutiveLaw::Clone() const
{
    return make_shared<PiecewiseLinearMomentCapacityConstitutiveLaw>(*this);
}

void PiecewiseLinearMomentCapacityConstitutiveLaw::GetLawFeatures(Features& rFeatures)
{
    rFeatures.mStrainSize     = strain_size;
    rFeatures.mSpaceDimension = space_dimenstion;
}

double& PiecewiseLinearMomentCapacityConstitutiveLaw::CalculateValue(ConstitutiveLaw::Parameters& rParameters,
                                                                     const Variable<double>& rVariable,
                                                                     double& rValue)
{
    if (rVariable == TANGENT_MODULUS) {
        const auto abs_curvature = std::abs(rParameters.GetStrainVector()[0]);
        // TANGENT_MODULUS is the derivative dM/d(kappa)
        rValue = mStressStrainTable.GetDerivative(abs_curvature);
        return rValue;
    }

    return BaseType::CalculateValue(rParameters, rVariable, rValue);
}

void PiecewiseLinearMomentCapacityConstitutiveLaw::CalculateMaterialResponsePK2(Parameters& rValues)
{
    const auto curvature = rValues.GetStrainVector()[0];

    const double moment          = CalculateMomentResponse(std::abs(curvature));
    auto&        r_stress_vector = rValues.GetStressVector();
    r_stress_vector[0]           = curvature < 0.0 ? -moment : moment;
}

SizeType PiecewiseLinearMomentCapacityConstitutiveLaw::GetStrainSize() const { return strain_size; }

int PiecewiseLinearMomentCapacityConstitutiveLaw::Check(const Properties&   rMaterialProperties,
                                                        const GeometryType& rElementGeometry,
                                                        const ProcessInfo& rCurrentProcessInfo) const
{
    if (const auto exit_code = BaseType::Check(rMaterialProperties, rElementGeometry, rCurrentProcessInfo);
        exit_code != 0)
        return exit_code;

    KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(STRAINS_OF_PIECEWISE_LINEAR_LAW))
        << "No STRAINS_OF_PIECEWISE_LINEAR_LAW found" << std::endl;
    KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(STRESSES_OF_PIECEWISE_LINEAR_LAW))
        << "No STRESSES_OF_PIECEWISE_LINEAR_LAW found" << std::endl;

    KRATOS_ERROR_IF_NOT(rMaterialProperties[STRAINS_OF_PIECEWISE_LINEAR_LAW].size() ==
                        rMaterialProperties[STRESSES_OF_PIECEWISE_LINEAR_LAW].size())
        << "The number of strain components does not match the number of stress components" << std::endl;
    KRATOS_ERROR_IF(rMaterialProperties[STRAINS_OF_PIECEWISE_LINEAR_LAW].empty())
        << "STRAINS_OF_PIECEWISE_LINEAR_LAW is empty";

    // The diagram implicitly starts at (0,0); the user must provide at least one point
    KRATOS_ERROR_IF(rMaterialProperties[STRAINS_OF_PIECEWISE_LINEAR_LAW].size() < 1)
        << "STRAINS_OF_PIECEWISE_LINEAR_LAW must contain at least one point";

    // Since (0,0) is assumed, the first provided point must be non-zero
    const double     first_k = rMaterialProperties[STRAINS_OF_PIECEWISE_LINEAR_LAW][0];
    const double     first_m = rMaterialProperties[STRESSES_OF_PIECEWISE_LINEAR_LAW][0];
    constexpr double tol     = 1.0e-15;
    KRATOS_ERROR_IF(std::abs(first_k) < tol || std::abs(first_m) < tol)
        << "First provided point must be non-zero when assuming implicit (0,0). Got (" << first_k
        << ", " << first_m << ")" << std::endl;

    // Allow non-decreasing sequence for piecewise M-k diagram (horizontal segments allowed)
    constexpr auto allow_equal = true;
    CheckUtilities::CheckValuesAreAscending(rMaterialProperties[STRAINS_OF_PIECEWISE_LINEAR_LAW],
                                            "STRAINS_OF_PIECEWISE_LINEAR_LAW", allow_equal);

    // Ensure the provided moment capacities are also non-decreasing with curvature
    CheckUtilities::CheckValuesAreAscending(rMaterialProperties[STRESSES_OF_PIECEWISE_LINEAR_LAW],
                                            "STRESSES_OF_PIECEWISE_LINEAR_LAW", allow_equal);

    return 0;
}

std::string PiecewiseLinearMomentCapacityConstitutiveLaw::Info() const
{
    return "PiecewiseLinearMomentCapacityConstitutiveLaw"s;
}

double PiecewiseLinearMomentCapacityConstitutiveLaw::CalculateMomentResponse(double AbsCurvature) const
{
    return mStressStrainTable.GetValue(AbsCurvature);
}

void PiecewiseLinearMomentCapacityConstitutiveLaw::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseType)
    rSerializer.save("StressStrainTable", mStressStrainTable);
}

void PiecewiseLinearMomentCapacityConstitutiveLaw::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseType)
    rSerializer.load("StressStrainTable", mStressStrainTable);
}

void PiecewiseLinearMomentCapacityConstitutiveLaw::InitializeMaterial(const Properties& rMaterialProperties,
                                                                      const ConstitutiveLaw::GeometryType& rElementGeometry,
                                                                      const Vector& rShapeFunctionsValues)
{
    BaseType::InitializeMaterial(rMaterialProperties, rElementGeometry, rShapeFunctionsValues);

    mStressStrainTable.Clear();

    const auto& r_strains = rMaterialProperties[STRAINS_OF_PIECEWISE_LINEAR_LAW];
    const auto& r_moments = rMaterialProperties[STRESSES_OF_PIECEWISE_LINEAR_LAW];
    // include implicit origin (0,0) as the first table point
    mStressStrainTable.PushBack(0.0, 0.0);
    for (auto i = std::size_t{0}; i < r_strains.size(); ++i) {
        mStressStrainTable.PushBack(r_strains[i], r_moments[i]);
    }
    // Append a final plateau point beyond the last provided strain so the
    // piecewise-linear table remains constant for larger curvatures
    if (!r_strains.empty()) {
        const std::size_t last_index = r_strains.size() - 1;
        const double      last_k     = r_strains[last_index];
        const double      last_m     = r_moments[last_index];
        mStressStrainTable.PushBack(last_k + 1.0, last_m);
    }
}

} // namespace Kratos
