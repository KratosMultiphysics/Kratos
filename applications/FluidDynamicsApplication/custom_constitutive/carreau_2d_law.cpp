//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Nicolò Antonelli
//

#include <iostream>
#include <cmath>
#include "includes/properties.h"
#include "custom_constitutive/carreau_2d_law.h"
#include "fluid_dynamics_application_variables.h"
#include "includes/cfd_variables.h"
#include "includes/checks.h"

namespace Kratos
{

Carreau2DLaw::Carreau2DLaw()
    : FluidConstitutiveLaw()
{
}

Carreau2DLaw::Carreau2DLaw(const Carreau2DLaw& rOther)
    : FluidConstitutiveLaw(rOther)
{
}

ConstitutiveLaw::Pointer Carreau2DLaw::Clone() const
{
    return Kratos::make_shared<Carreau2DLaw>(*this);
}

Carreau2DLaw::~Carreau2DLaw()
{
}

ConstitutiveLaw::SizeType Carreau2DLaw::WorkingSpaceDimension()
{
    return 2;
}

ConstitutiveLaw::SizeType Carreau2DLaw::GetStrainSize() const
{
    return 3;
}

void Carreau2DLaw::CalculateMaterialResponseCauchy(Parameters& rValues)
{
    Flags& Options = rValues.GetOptions();
    const Properties& MaterialProperties = rValues.GetMaterialProperties();

    Vector& S = rValues.GetStrainVector(); // Strain rate tensor
    Vector& StressVector = rValues.GetStressVector();

    // Material parameters
    const double eta_0 = MaterialProperties[ZERO_SHEAR_VISCOSITY];
    const double eta_inf = MaterialProperties[INFINITE_SHEAR_VISCOSITY];
    const double lambda = MaterialProperties[RELAXATION_TIME];
    const double n = MaterialProperties[SHEAR_THINNING_INDEX];
    const double a = MaterialProperties[CARREAU_TRANSITION_SHARPNESS];

    // Compute gamma_dot (strain rate magnitude)
    const double gamma_dot = std::sqrt(2.0 * S[0] * S[0] + 2.0 * S[1] * S[1] + S[2] * S[2]);
    const double min_gamma_dot = 1e-12; // Regularization for zero strain rate
    const double g = std::max(gamma_dot, min_gamma_dot);

    // Compute effective viscosity based on Carreau model
    const double eta_effective = eta_inf + (eta_0 - eta_inf) * std::pow(1.0 + std::pow(lambda * g, a), (n - 1.0) / a);
    // KRATOS_WATCH(eta_effective)
    // Compute stress tensor
    StressVector[0] = 2.0 * eta_effective * S[0];
    StressVector[1] = 2.0 * eta_effective * S[1];
    StressVector[2] = eta_effective * S[2];

    if (Options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR))
    {
        this->NewtonianConstitutiveMatrix2D(eta_effective, rValues.GetConstitutiveMatrix());
    }
}

int Carreau2DLaw::Check(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const ProcessInfo& rCurrentProcessInfo) const
{
    if (rMaterialProperties[ZERO_SHEAR_VISCOSITY] <= 0.0)
    {
        KRATOS_ERROR << "ZERO_SHEAR_VISCOSITY must be greater than zero." << std::endl;
    }

    if (rMaterialProperties[INFINITE_SHEAR_VISCOSITY] < 0.0)
    {
        KRATOS_ERROR << "INFINITE_SHEAR_VISCOSITY cannot be negative." << std::endl;
    }

    if (rMaterialProperties[RELAXATION_TIME] <= 0.0)
    {
        KRATOS_ERROR << "RELAXATION_TIME must be greater than zero." << std::endl;
    }

    if (rMaterialProperties[SHEAR_THINNING_INDEX] <= 0.0 || rMaterialProperties[SHEAR_THINNING_INDEX] > 1.0)
    {
        KRATOS_ERROR << "SHEAR_THINNING_INDEX must be in the range (0, 1]." << std::endl;
    }

    if (rMaterialProperties[CARREAU_TRANSITION_SHARPNESS] <= 0.0)
    {
        KRATOS_ERROR << "CARREAU_TRANSITION_SHARPNESS must be greater than zero." << std::endl;
    }

    return 0;
}

double Carreau2DLaw::GetEffectiveViscosity(ConstitutiveLaw::Parameters& rParameters) const {
    // We are abusing the fact that C(2,2) = mu_effective
    return rParameters.GetConstitutiveMatrix()(2,2);
}

std::string Carreau2DLaw::Info() const
{
    return "Carreau2DLaw";
}

void Carreau2DLaw::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, FluidConstitutiveLaw)
}

void Carreau2DLaw::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, FluidConstitutiveLaw)
}

} // namespace Kratos
