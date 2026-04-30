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
        const auto moment_response =
            CalculateAbsMomentResponse(abs_curvature, rParameters.GetMaterialProperties());
        rValue = moment_response.tangent_modulus;
    } else {
        rValue = BaseType::CalculateValue(rParameters, rVariable, rValue);
    }

    return rValue;
}

void PiecewiseLinearMomentCapacityConstitutiveLaw::CalculateMaterialResponsePK2(Parameters& rValues)
{
    const auto curvature = rValues.GetStrainVector()[0];

    const auto moment_response =
        CalculateAbsMomentResponse(std::abs(curvature), rValues.GetMaterialProperties());

    auto& r_stress_vector = rValues.GetStressVector();
    r_stress_vector[0]    = curvature < 0.0 ? -moment_response.moment : moment_response.moment;
}

SizeType PiecewiseLinearMomentCapacityConstitutiveLaw::GetStrainSize() const { return strain_size; }

int PiecewiseLinearMomentCapacityConstitutiveLaw::Check(const Properties&   rMaterialProperties,
                                                        const GeometryType& rElementGeometry,
                                                        const ProcessInfo& rCurrentProcessInfo) const
{
    if (const auto exit_code = BaseType::Check(rMaterialProperties, rElementGeometry, rCurrentProcessInfo);
        exit_code != 0)
        return exit_code;

    const auto checker =
        CheckProperties{rMaterialProperties, "Properties", CheckProperties::Bounds::AllInclusive};
    checker.CheckAvailability(STRAINS_OF_PIECEWISE_LINEAR_LAW);
    checker.CheckAvailability(STRESSES_OF_PIECEWISE_LINEAR_LAW);
    checker.Check(MAX_AXIAL_LOAD_OF_CONSTRUCTION_ELEMENT);
    checker.Check(MOMENT_CAPACITY_REDUCTION_FACTOR);
    checker.Check(MINIMUM_MOMENT_CAPACITY_FACTOR, 1.0);

    const auto& r_curvatures        = rMaterialProperties[STRAINS_OF_PIECEWISE_LINEAR_LAW];
    const auto& r_moment_capacities = rMaterialProperties[STRESSES_OF_PIECEWISE_LINEAR_LAW];

    // Array size constants
    static constexpr SizeType number_of_curvature_points  = 3;
    static constexpr SizeType number_of_moment_capacities = 2;

    KRATOS_ERROR_IF_NOT(r_curvatures.size() == number_of_curvature_points)
        << "STRAINS_OF_PIECEWISE_LINEAR_LAW requires exactly " << number_of_curvature_points
        << " values: [kappa_elastic_end, kappa_plastic1_end, "
           "kappa_elastoplastic_end]"
        << std::endl;
    KRATOS_ERROR_IF_NOT(r_moment_capacities.size() == number_of_moment_capacities)
        << "STRESSES_OF_PIECEWISE_LINEAR_LAW requires exactly " << number_of_moment_capacities
        << " values: [M_plastic1, M_plateau2]" << std::endl;

    CheckCurvaturesAreAscending(r_curvatures);

    KRATOS_ERROR_IF(r_curvatures[0] <= 0.0) << "STRAINS_OF_PIECEWISE_LINEAR_LAW[0] must be > 0.0" << std::endl;
    KRATOS_ERROR_IF(r_moment_capacities[0] <= 0.0)
        << "STRESSES_OF_PIECEWISE_LINEAR_LAW[0] must be > 0.0" << std::endl;
    KRATOS_ERROR_IF(r_moment_capacities[1] < r_moment_capacities[0])
        << "STRESSES_OF_PIECEWISE_LINEAR_LAW[1] must be >= STRESSES_OF_PIECEWISE_LINEAR_LAW[0]"
        << std::endl;

    return 0;
}

std::string PiecewiseLinearMomentCapacityConstitutiveLaw::Info() const
{
    return "PiecewiseLinearMomentCapacityConstitutiveLaw"s;
}

std::pair<double, double> PiecewiseLinearMomentCapacityConstitutiveLaw::CalculateScaledMomentCapacities(const Properties& rMaterialProperties)
{
    const auto& r_moment_capacities = rMaterialProperties[STRESSES_OF_PIECEWISE_LINEAR_LAW];

    const auto max_axial_load = rMaterialProperties[MAX_AXIAL_LOAD_OF_CONSTRUCTION_ELEMENT];
    const auto reduction      = rMaterialProperties[MOMENT_CAPACITY_REDUCTION_FACTOR];
    const auto min_factor     = rMaterialProperties[MINIMUM_MOMENT_CAPACITY_FACTOR];

    const auto reduction_factor = std::clamp(1.0 - reduction * max_axial_load, min_factor, 1.0);

    return {reduction_factor * r_moment_capacities[0], reduction_factor * r_moment_capacities[1]};
}

MomentResponse PiecewiseLinearMomentCapacityConstitutiveLaw::CalculateAbsMomentResponse(double AbsCurvature,
                                                                                        const Properties& rMaterialProperties)
{
    const auto& r_curvatures        = rMaterialProperties[STRAINS_OF_PIECEWISE_LINEAR_LAW];
    const auto [moment_1, moment_2] = CalculateScaledMomentCapacities(rMaterialProperties);

    const auto curvature_1 = r_curvatures[0];
    const auto curvature_2 = r_curvatures[1];
    const auto curvature_3 = r_curvatures[2];

    if (AbsCurvature <= curvature_1) {
        const auto initial_slope = moment_1 / curvature_1;
        return {initial_slope * AbsCurvature, initial_slope};
    }

    if (AbsCurvature <= curvature_2) {
        return {moment_1, 0.0};
    }

    if (AbsCurvature <= curvature_3) {
        const auto hardening_slope = (moment_2 - moment_1) / (curvature_3 - curvature_2);
        return {moment_1 + hardening_slope * (AbsCurvature - curvature_2), hardening_slope};
    }

    return {moment_2, 0.0};
}

void PiecewiseLinearMomentCapacityConstitutiveLaw::CheckCurvaturesAreAscending(const Vector& rCurvatures)
{
    auto first_ge_second = [](const auto& First, const auto& Second) { return First >= Second; };
    auto pos = std::adjacent_find(rCurvatures.cbegin(), rCurvatures.cend(), first_ge_second);
    KRATOS_ERROR_IF(pos != rCurvatures.cend())
        << "Values in STRAINS_OF_PIECEWISE_LINEAR_LAW are not ascending: " << *(pos + 1)
        << " (at index " << std::distance(rCurvatures.begin(), pos) + 2 << ") does not exceed "
        << *pos << " (at index " << std::distance(rCurvatures.begin(), pos) + 1 << ")" << std::endl;
}

void PiecewiseLinearMomentCapacityConstitutiveLaw::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseType)
}

void PiecewiseLinearMomentCapacityConstitutiveLaw::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseType)
}

} // namespace Kratos
