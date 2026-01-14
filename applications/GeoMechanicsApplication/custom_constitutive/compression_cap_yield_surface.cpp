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
#include "custom_utilities/check_utilities.h"
#include "custom_utilities/constitutive_law_utilities.h"
#include "custom_utilities/ublas_utilities.h"
#include "geo_mechanics_application_variables.h"
#include "includes/serializer.h"
#include "utilities/string_utilities.h"

#include <cmath>

namespace
{

using namespace Kratos;

double GetCapSize(const Properties& rProperties)
{
    if (rProperties.Has(GEO_COMPRESSION_CAP_SIZE)) {
        return rProperties[GEO_COMPRESSION_CAP_SIZE];
    }
    if (rProperties.Has(K0_NC)) {
        return 4.5 * (1.0 - rProperties[K0_NC]) / (1.0 + 2.0 * rProperties[K0_NC]);
    }
    if (rProperties.Has(GEO_FRICTION_ANGLE)) {
        auto k0_nc = 1.0 - std::sin(ConstitutiveLawUtilities::GetFrictionAngleInRadians(rProperties));
        return 4.5 * (1.0 - k0_nc) / (1.0 + 2.0 * k0_nc);
    }
    KRATOS_ERROR << "ConstitutiveLawUtilities::GetCapSize failed. There is no "
                    "GEO_COMPRESSION_CAP_SIZE available "
                 << std::endl;
}

double GetCapLocation(const Properties& rProperties)
{
    if (rProperties.Has(GEO_COMPRESSION_CAP_LOCATION)) {
        return rProperties[GEO_COMPRESSION_CAP_LOCATION];
    }
    KRATOS_ERROR << "ConstitutiveLawUtilities::GetCapLocation failed. There is no "
                    "GEO_COMPRESSION_CAP_LOCATION available "
                 << std::endl;
}

CompressionCapYieldSurface::KappaDependentFunction MakeConstantCapFunction(double Value)
{
    return [Value](double /* unused kappa */) { return Value; };
}

std::string GetCapHardeningTypeFrom(const Properties& rMaterialProperties)
{
    return StringUtilities::ConvertCamelCaseToSnakeCase(rMaterialProperties[GEO_CAP_HARDENING_TYPE]);
}

CompressionCapYieldSurface::KappaDependentFunction MakeCapSizeCalculator(const Properties& rMaterialProperties)
{
    const auto hardening_type = GetCapHardeningTypeFrom(rMaterialProperties);
    if (hardening_type == "none") {
        return MakeConstantCapFunction(GetCapSize(rMaterialProperties));
    }
    KRATOS_ERROR << "Cannot create a kappa-dependent function for the cap size of material "
                 << rMaterialProperties.Id() << ": unknown hardening type '" << hardening_type << "'\n";
}

CompressionCapYieldSurface::KappaDependentFunction MakeCapLocationCalculator(const Properties& rMaterialProperties)
{
    const auto hardening_type = GetCapHardeningTypeFrom(rMaterialProperties);
    if (hardening_type == "none") {
        return MakeConstantCapFunction(GetCapLocation(rMaterialProperties));
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

double CompressionCapYieldSurface::YieldFunctionValue(const Vector& rSigmaTau) const
{
    return std::pow(rSigmaTau[1] / GetCapSize(), 2) +
           (rSigmaTau[0] + GetCapLocation()) * (rSigmaTau[0] - GetCapLocation());
}

Vector CompressionCapYieldSurface::DerivativeOfFlowFunction(const Vector& rSigmaTau) const
{
    return UblasUtilities::CreateVector({2.0 * rSigmaTau[0], 2.0 * rSigmaTau[1] / std::pow(GetCapSize(), 2)});
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
