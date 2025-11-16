//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

// System includes
#include <iostream>

// External includes

// Project includes
#include "includes/cfd_variables.h"
#include "includes/checks.h"

// Application includes
#include "custom_constitutive/newtonian_2d_law.h"
#include "custom_constitutive/newtonian_3d_law.h"
#include "custom_utilities/fluid_calculation_utilities.h"
#include "rans_application_variables.h"

// Include base h
#include "custom_constitutive/rans_k_epsilon_newtonian_law.h"

namespace Kratos
{
//******************************CONSTRUCTOR*******************************************
//************************************************************************************

template<class TPrimalBaseType>
RansKEpsilonNewtonianLaw<TPrimalBaseType>::RansKEpsilonNewtonianLaw() : BaseType()
{
}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

template<class TPrimalBaseType>
RansKEpsilonNewtonianLaw<TPrimalBaseType>::RansKEpsilonNewtonianLaw(const RansKEpsilonNewtonianLaw& rOther)
    : BaseType(rOther)
{
}

//********************************CLONE***********************************************
//************************************************************************************

template<class TPrimalBaseType>
ConstitutiveLaw::Pointer RansKEpsilonNewtonianLaw<TPrimalBaseType>::Clone() const
{
    return Kratos::make_shared<RansKEpsilonNewtonianLaw>(*this);
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

template<class TPrimalBaseType>
RansKEpsilonNewtonianLaw<TPrimalBaseType>::~RansKEpsilonNewtonianLaw()
{
}

template<class TPrimalBaseType>
int RansKEpsilonNewtonianLaw<TPrimalBaseType>::Check(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    // Check viscosity value
    KRATOS_ERROR_IF(rMaterialProperties[DYNAMIC_VISCOSITY] <= 0.0)
        << "Incorrect or missing DYNAMIC_VISCOSITY provided in material properties "
           "for RansKEpsilonNewtonianLaw: "
        << rMaterialProperties[DYNAMIC_VISCOSITY] << std::endl;

    KRATOS_ERROR_IF(rMaterialProperties[DENSITY] <= 0.0)
        << "Incorrect or missing DENSITY provided in material properties "
           "for RansKEpsilonNewtonianLaw: "
        << rMaterialProperties[DENSITY] << std::endl;

    KRATOS_ERROR_IF(rCurrentProcessInfo[TURBULENCE_RANS_C_MU] <= 0.0)
        << "Incorrect or missing TURBULENCE_RANS_C_MU provided in process info "
           "for RansKEpsilonNewtonianLaw: "
        << rCurrentProcessInfo[TURBULENCE_RANS_C_MU] << std::endl;

    return 0;

    KRATOS_CATCH("");
}

template<class TPrimalBaseType>
double& RansKEpsilonNewtonianLaw<TPrimalBaseType>::CalculateValue(
    ConstitutiveLaw::Parameters& rParameters,
    const Variable<double>& rVariable,
    double& rValue)
{
    KRATOS_TRY

    if (rVariable == TURBULENT_VISCOSITY) {
        rValue = this->CalculateTurbulentViscosity(rParameters);
    } else if (rVariable == DYNAMIC_VISCOSITY) {
        rValue = this->GetDynamicViscosity(rParameters);
    } else {
        rValue = BaseType::CalculateValue(rParameters, rVariable, rValue);
    }

    return rValue;

    KRATOS_CATCH("");
}

template<class TPrimalBaseType>
void RansKEpsilonNewtonianLaw<TPrimalBaseType>::CalculateDerivative(
    ConstitutiveLaw::Parameters& rParameters,
    const Variable<double>& rFunctionVariable,
    const Variable<double>& rDerivativeVariable,
    double& rOutput)
{
    KRATOS_TRY

    if (rFunctionVariable == EFFECTIVE_VISCOSITY) {
        const double c_mu = rParameters.GetProcessInfo()[TURBULENCE_RANS_C_MU];
        const Properties& r_prop = rParameters.GetMaterialProperties();
        const double density = r_prop[DENSITY];

        rOutput = 0.0;

        double tke, epsilon;
        FluidCalculationUtilities::EvaluateInPoint(
            rParameters.GetElementGeometry(), rParameters.GetShapeFunctionsValues(),
            std::tie(tke, TURBULENT_KINETIC_ENERGY),
            std::tie(epsilon, TURBULENT_ENERGY_DISSIPATION_RATE));

        const double nu_t = c_mu * std::pow(tke, 2) / epsilon;
        if (epsilon > 0.0 && nu_t > 1e-12) {
            // computing gauss point nu_t derivative
            if (rDerivativeVariable == TURBULENT_KINETIC_ENERGY) {
                rOutput = 2.0 * c_mu * tke / epsilon;
            } else if (rDerivativeVariable == TURBULENT_ENERGY_DISSIPATION_RATE) {
                rOutput = -1.0 * c_mu * std::pow(tke / epsilon, 2);
            }
        }

        rOutput *= density;
    }

    KRATOS_CATCH("");
}

template<class TPrimalBaseType>
std::string RansKEpsilonNewtonianLaw<TPrimalBaseType>::Info() const
{
    return "RansKEpsilon" + TPrimalBaseType::Info();
}

template<class TPrimalBaseType>
double RansKEpsilonNewtonianLaw<TPrimalBaseType>::GetEffectiveViscosity(ConstitutiveLaw::Parameters& rParameters) const
{
    const Properties& r_prop = rParameters.GetMaterialProperties();
    const double mu = r_prop[DYNAMIC_VISCOSITY];
    const double density = r_prop[DENSITY];

    return mu + density * this->CalculateTurbulentViscosity(rParameters);
}

template<class TPrimalBaseType>
double RansKEpsilonNewtonianLaw<TPrimalBaseType>::CalculateTurbulentViscosity(
    ConstitutiveLaw::Parameters& rParameters) const
{
    const double c_mu = rParameters.GetProcessInfo()[TURBULENCE_RANS_C_MU];

    double tke, epsilon;
    FluidCalculationUtilities::EvaluateInPoint(
        rParameters.GetElementGeometry(), rParameters.GetShapeFunctionsValues(),
        std::tie(tke, TURBULENT_KINETIC_ENERGY),
        std::tie(epsilon, TURBULENT_ENERGY_DISSIPATION_RATE));

    double nu_t = 1e-12;
    if (epsilon > 0.0) {
        nu_t = std::max(c_mu * std::pow(tke, 2) / epsilon, 1e-12);
    }

    return nu_t;
}

template<class TPrimalBaseType>
double RansKEpsilonNewtonianLaw<TPrimalBaseType>::GetDynamicViscosity(
    ConstitutiveLaw::Parameters& rParameters) const
{
    const Properties& r_prop = rParameters.GetMaterialProperties();
    return r_prop[DYNAMIC_VISCOSITY];
}

template<class TPrimalBaseType>
void RansKEpsilonNewtonianLaw<TPrimalBaseType>::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseType)
}

template<class TPrimalBaseType>
void RansKEpsilonNewtonianLaw<TPrimalBaseType>::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseType)
}

// template instantiations
template class RansKEpsilonNewtonianLaw<Newtonian2DLaw>;
template class RansKEpsilonNewtonianLaw<Newtonian3DLaw>;

} // Namespace Kratos
