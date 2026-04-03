// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Mohamed Nabi
//                   Wijtze Pieter Kikstra
//                   Anne van de Graaf
//

#include "custom_constitutive/compression_cap_yield_surface.h"
#include "custom_constitutive/p_q.hpp"
#include "custom_utilities/check_utilities.hpp"
#include "custom_utilities/constitutive_law_utilities.h"
#include "custom_utilities/function_object_utilities.h"
#include "custom_utilities/math_utilities.hpp"
#include "custom_utilities/stress_strain_utilities.h"
#include "custom_utilities/string_utilities.h"
#include "custom_utilities/ublas_utilities.h"
#include "geo_mechanics_application_variables.h"
#include "includes/serializer.h"

#include <cmath>

namespace
{

using namespace Kratos;

double CalculateCapSizeFromK0NC(double K0_nc) { return 4.5 * (1.0 - K0_nc) / (1.0 + 2.0 * K0_nc); }

double GetCapSize(const Properties& rProperties)
{
    if (rProperties.Has(GEO_COMPRESSION_CAP_SIZE)) {
        return rProperties[GEO_COMPRESSION_CAP_SIZE];
    }

    if (rProperties.Has(K0_NC)) {
        return CalculateCapSizeFromK0NC(rProperties[K0_NC]);
    }

    if (rProperties.Has(GEO_FRICTION_ANGLE)) {
        const auto k0_nc = ConstitutiveLawUtilities::CalculateK0NCFromFrictionAngleInRadians(
            MathUtils<>::DegreesToRadians(rProperties[GEO_FRICTION_ANGLE]));
        return CalculateCapSizeFromK0NC(k0_nc);
    }

    KRATOS_ERROR
        << "Failed to determine the compression cap size. Neither of the following material "
           "properties was provided: GEO_COMPRESSION_CAP_SIZE, K0_NC, or GEO_FRICTION_ANGLE\n ";
}

double GetPreconsolidationStress(const Properties& rProperties)
{
    if (rProperties.Has(GEO_PRECONSOLIDATION_STRESS)) {
        return rProperties[GEO_PRECONSOLIDATION_STRESS];
    }
    KRATOS_ERROR << "Failed to determine the preconsolidation stress. There is no "
                    "GEO_PRECONSOLIDATION_STRESS available "
                 << std::endl;
}

Geo::KappaDependentFunction MakeCapSizeCalculator(const Properties& rMaterialProperties)
{
    const auto hardening_type = GeoStringUtilities::ToLower(rMaterialProperties[GEO_CAP_HARDENING_TYPE]);
    if (hardening_type == "none") {
        return FunctionObjectUtilities::MakeConstantFunction(GetCapSize(rMaterialProperties));
    }
    KRATOS_ERROR << "Cannot create a kappa-dependent function for the cap size of material "
                 << rMaterialProperties.Id() << ": unknown hardening type '" << hardening_type << "'\n";
}

Geo::KappaDependentFunction MakeCapLocationCalculator(const Properties& rMaterialProperties)
{
    const auto hardening_type = GeoStringUtilities::ToLower(rMaterialProperties[GEO_CAP_HARDENING_TYPE]);
    if (hardening_type == "none") {
        return FunctionObjectUtilities::MakeConstantFunction(GetPreconsolidationStress(rMaterialProperties));
    }
    KRATOS_ERROR << "Cannot create a kappa-dependent function for the cap location of material "
                 << rMaterialProperties.Id() << ": unknown hardening type '" << hardening_type << "'\n";
}

} // namespace

namespace Kratos
{

CompressionCapYieldSurface::CompressionCapYieldSurface()
{
    mMaterialProperties.SetValue(GEO_CAP_HARDENING_TYPE, "None");
    mMaterialProperties.SetValue(GEO_COMPRESSION_CAP_SIZE, 0.0);
    mMaterialProperties.SetValue(GEO_PRECONSOLIDATION_STRESS, 0.0);

    InitializeKappaDependentFunctions();
}

CompressionCapYieldSurface::CompressionCapYieldSurface(const Properties& rMaterialProperties)
    : mMaterialProperties{rMaterialProperties}
{
    if (!mMaterialProperties.Has(GEO_CAP_HARDENING_TYPE)) {
        mMaterialProperties.SetValue(GEO_CAP_HARDENING_TYPE, "None");
    }

    CheckMaterialProperties();
    InitializeKappaDependentFunctions();
}

double CompressionCapYieldSurface::GetCapSize() const { return mCapSizeCalculator(mKappa); }

double CompressionCapYieldSurface::GetPreconsolidationStress() const
{
    return mPreconsolidationStressCalculator(mKappa);
}

double CompressionCapYieldSurface::YieldFunctionValue(const Geo::PQ& rPQ) const
{
    return std::pow(rPQ.Q() / GetCapSize(), 2) + rPQ.P() * rPQ.P() - std::pow(GetPreconsolidationStress(), 2);
}

double CompressionCapYieldSurface::YieldFunctionValue(const Geo::PrincipalStresses& rPrincipalStresses) const
{
    return YieldFunctionValue(StressStrainUtilities::TransformPrincipalStressesToPandQ(rPrincipalStresses));
}

double CompressionCapYieldSurface::YieldFunctionValue(const Geo::SigmaTau&) const
{
    KRATOS_ERROR << "CompressionCapYieldSurface::YieldFunctionValue(const Geo::SigmaTau&) is not "
                    "implemented\n";
}

Vector CompressionCapYieldSurface::DerivativeOfFlowFunction(const Geo::PQ& rPQ) const
{
    return UblasUtilities::CreateVector({2.0 * rPQ.P(), 2.0 * rPQ.Q() / std::pow(GetCapSize(), 2)});
}

Vector CompressionCapYieldSurface::DerivativeOfFlowFunction(const Geo::PrincipalStresses& rPrincipalStresses) const
{
    const auto p_q = StressStrainUtilities::TransformPrincipalStressesToPandQ(rPrincipalStresses);
    auto       derivative_pq = DerivativeOfFlowFunction(p_q);
    const auto c1            = derivative_pq[0] / 3.0;
    const auto c2            = 3.0 * derivative_pq[1] / (2.0 * p_q.Q());
    Vector     result        = Vector(3, 0.0);
    result[0]                = c1 + c2 * (rPrincipalStresses.Values()[0] - p_q.P());
    result[1]                = c1 + c2 * (rPrincipalStresses.Values()[1] - p_q.P());
    result[2]                = c1 + c2 * (rPrincipalStresses.Values()[2] - p_q.P());
    return result;
}

void CompressionCapYieldSurface::InitializeKappaDependentFunctions()
{
    mCapSizeCalculator                = MakeCapSizeCalculator(mMaterialProperties);
    mPreconsolidationStressCalculator = MakeCapLocationCalculator(mMaterialProperties);
}

void CompressionCapYieldSurface::CheckMaterialProperties() const
{
    const CheckProperties check_properties(mMaterialProperties, "property", CheckProperties::Bounds::AllInclusive);
    check_properties.Check(GEO_PRECONSOLIDATION_STRESS);
    if (mMaterialProperties.Has(GEO_COMPRESSION_CAP_SIZE)) {
        check_properties.Check(GEO_COMPRESSION_CAP_SIZE);
    } else if (mMaterialProperties.Has(K0_NC)) {
        check_properties.Check(K0_NC);
    } else if (mMaterialProperties.Has(GEO_FRICTION_ANGLE)) {
        check_properties.Check(GEO_FRICTION_ANGLE);
    } else {
        KRATOS_ERROR << "Invalid material properties: cap size determination requires one of "
                        "GEO_COMPRESSION_CAP_SIZE, K0_NC, or GEO_FRICTION_ANGLE to be defined."
                     << std::endl;
    }
}

double CompressionCapYieldSurface::CalculatePlasticMultiplier(const Geo::PrincipalStresses& rPrincipalStresses,
                                                              const Vector& rDerivativeOfFlowFunction,
                                                              const Matrix& rElasticConstitutiveTensor) const
{
    const auto principals =
        UblasUtilities::CreateVector({rPrincipalStresses.Values()[0] - rPrincipalStresses.Values()[1],
                                      rPrincipalStresses.Values()[1] - rPrincipalStresses.Values()[2],
                                      rPrincipalStresses.Values()[2] - rPrincipalStresses.Values()[0]});
    const auto elastic_matrix    = subrange(rElasticConstitutiveTensor, 0, 3, 0, 3);
    const auto stress_correction = Vector{prod(elastic_matrix, rDerivativeOfFlowFunction)};
    const auto temp = UblasUtilities::CreateVector({stress_correction[0] - stress_correction[1],
                                                    stress_correction[1] - stress_correction[2],
                                                    stress_correction[2] - stress_correction[0]});
    const auto A    = inner_prod(temp, temp) / (2.0 * std::pow(GetCapSize(), 2)) +
                   std::pow((stress_correction[0] + stress_correction[1] + stress_correction[2]) / 3.0, 2);
    const auto B =
        inner_prod(principals, temp) / std::pow(GetCapSize(), 2) +
        (rPrincipalStresses.Values()[0] + rPrincipalStresses.Values()[1] + rPrincipalStresses.Values()[2]) *
            (stress_correction[0] + stress_correction[1] + stress_correction[2]) * 2.0 / 9.0;
    const auto C = YieldFunctionValue(rPrincipalStresses);

    const auto roots = GeoMechanicsMathUtilities::RootsOfSecondOrderEquation(A, B, C);
    KRATOS_DEBUG_ERROR_IF(roots.size() == 0)
        << "Failed to calculate the plastic multiplier for the cap.\n";
    return *std::ranges::max_element(roots);
}

void CompressionCapYieldSurface::save(Serializer& rSerializer) const
{
    rSerializer.save("Kappa", mKappa);
    rSerializer.save("MaterialProperties", mMaterialProperties);
}

void CompressionCapYieldSurface::load(Serializer& rSerializer)
{
    rSerializer.load("Kappa", mKappa);
    rSerializer.load("MaterialProperties", mMaterialProperties);

    InitializeKappaDependentFunctions();
}

} // Namespace Kratos
