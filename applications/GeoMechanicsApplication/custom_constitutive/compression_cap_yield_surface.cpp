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
//

#include "custom_constitutive/compression_cap_yield_surface.h"
#include "custom_utilities/check_utilities.h"
#include "custom_utilities/constitutive_law_utilities.h"
#include "custom_utilities/ublas_utilities.h"
#include "geo_mechanics_application_variables.h"
#include "includes/serializer.h"

#include <boost/numeric/ublas/assignment.hpp>
#include <cmath>

namespace
{

using namespace Kratos;

CompressionCapYieldSurface::KappaDependentFunction MakeConstantFunction(double Value)
{
    return [Value](double /* unused kappa */) { return Value; };
}

std::string GetCapHardeningTypeFrom(const Properties& rMaterialProperties)
{
    auto result   = rMaterialProperties[GEO_CAP_HARDENING_TYPE];
    auto to_lower = [](auto character) { return std::tolower(character); };
    std::ranges::transform(result, result.begin(), to_lower);
    return result;
}

CompressionCapYieldSurface::KappaDependentFunction MakeCapSizeCalculator(const Properties& rMaterialProperties)
{
    const auto hardening_type = GetCapHardeningTypeFrom(rMaterialProperties);
    if (hardening_type == "none") {
        return MakeConstantFunction(ConstitutiveLawUtilities::GetCapSize(rMaterialProperties));
    }
    KRATOS_ERROR << "Cannot create a kappa-dependent function for the cap size of material "
                 << rMaterialProperties.Id() << ": unknown hardening type '" << hardening_type << "'\n";
}

CompressionCapYieldSurface::KappaDependentFunction MakeCapLocationCalculator(const Properties& rMaterialProperties)
{
    const auto hardening_type = GetCapHardeningTypeFrom(rMaterialProperties);
    if (hardening_type == "none") {
        return MakeConstantFunction(ConstitutiveLawUtilities::GetCapLocation(rMaterialProperties));
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
    mMaterialProperties.SetValue(GEO_COMPRESSION_CAP_LOCATION, 0.0);

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

double CompressionCapYieldSurface::GetCapLocation() const { return mCapLocationCalculator(mKappa); }

double CompressionCapYieldSurface::GetKappa() const { return mKappa; }

void CompressionCapYieldSurface::SetKappa(double kappa) { mKappa = kappa; }

double CompressionCapYieldSurface::YieldFunctionValue(const Vector& rSigmaTau) const
{
    return std::pow(rSigmaTau[1] / GetCapSize(), 2) +
           (rSigmaTau[0] + GetCapLocation()) * (rSigmaTau[0] - GetCapLocation());
}

Vector CompressionCapYieldSurface::DerivativeOfFlowFunction(const Vector& rSigmaTau) const
{
    Vector result(2);
    result <<= 2.0 * rSigmaTau[0], 2.0 * rSigmaTau[1] / std::pow(GetCapSize(), 2);
    return result;
}

double CompressionCapYieldSurface::CalculateCapCornerPoint() const
{
    return 0; // TODO: implement this
}

void CompressionCapYieldSurface::InitializeKappaDependentFunctions()
{
    mCapSizeCalculator     = MakeCapSizeCalculator(mMaterialProperties);
    mCapLocationCalculator = MakeCapLocationCalculator(mMaterialProperties);
}

void CompressionCapYieldSurface::CheckMaterialProperties() const
{
    const CheckProperties check_properties(mMaterialProperties, "property", CheckProperties::Bounds::AllInclusive);
    check_properties.Check(GEO_COMPRESSION_CAP_SIZE);
    check_properties.Check(GEO_COMPRESSION_CAP_LOCATION);
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
