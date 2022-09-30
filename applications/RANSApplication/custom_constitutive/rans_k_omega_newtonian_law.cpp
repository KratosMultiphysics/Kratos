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
#include "custom_constitutive/rans_k_omega_newtonian_law.h"

namespace Kratos
{
//******************************CONSTRUCTOR*******************************************
//************************************************************************************

template<class TPrimalBaseType>
RansKOmegaNewtonianLaw<TPrimalBaseType>::RansKOmegaNewtonianLaw() : BaseType()
{
}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

template<class TPrimalBaseType>
RansKOmegaNewtonianLaw<TPrimalBaseType>::RansKOmegaNewtonianLaw(const RansKOmegaNewtonianLaw& rOther)
    : BaseType(rOther)
{
}

//********************************CLONE***********************************************
//************************************************************************************

template<class TPrimalBaseType>
ConstitutiveLaw::Pointer RansKOmegaNewtonianLaw<TPrimalBaseType>::Clone() const
{
    return Kratos::make_shared<RansKOmegaNewtonianLaw>(*this);
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

template<class TPrimalBaseType>
RansKOmegaNewtonianLaw<TPrimalBaseType>::~RansKOmegaNewtonianLaw()
{
}

template<class TPrimalBaseType>
int RansKOmegaNewtonianLaw<TPrimalBaseType>::Check(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    // Check viscosity value
    KRATOS_ERROR_IF(rMaterialProperties[DYNAMIC_VISCOSITY] <= 0.0)
        << "Incorrect or missing DYNAMIC_VISCOSITY provided in material properties "
           "for RansKOmegaNewtonianLaw: "
        << rMaterialProperties[DYNAMIC_VISCOSITY] << std::endl;

    KRATOS_ERROR_IF(rMaterialProperties[DENSITY] <= 0.0)
        << "Incorrect or missing DENSITY provided in material properties "
           "for RansKOmegaNewtonianLaw: "
        << rMaterialProperties[DENSITY] << std::endl;

    return 0;

    KRATOS_CATCH("");
}

template<class TPrimalBaseType>
double& RansKOmegaNewtonianLaw<TPrimalBaseType>::CalculateValue(
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
void RansKOmegaNewtonianLaw<TPrimalBaseType>::CalculateDerivative(
    ConstitutiveLaw::Parameters& rParameters,
    const Variable<double>& rFunctionVariable,
    const Variable<double>& rDerivativeVariable,
    double& rOutput)
{
    KRATOS_TRY

    if (rFunctionVariable == EFFECTIVE_VISCOSITY) {
        const Properties& r_prop = rParameters.GetMaterialProperties();
        const double density = r_prop[DENSITY];

        rOutput = 0.0;

        double tke, omega;
        FluidCalculationUtilities::EvaluateInPoint(
            rParameters.GetElementGeometry(), rParameters.GetShapeFunctionsValues(),
            std::tie(tke, TURBULENT_KINETIC_ENERGY),
            std::tie(omega, TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE));

        double nu_t = 1e-12;
        if (tke > 1e-12 && omega > 1e-12) {
            nu_t = std::max(tke / omega, 1e-12);
        }
        if (tke > 0.0 && omega > 1e-12 && nu_t > 1e-12) {
            // computing gauss point nu_t derivative
            if (rDerivativeVariable == TURBULENT_KINETIC_ENERGY) {
                rOutput = 1 / omega;
            } else if (rDerivativeVariable == TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE) {
                rOutput = -tke / std::pow(omega, 2);
            }
        }

        rOutput *= density;
    }

    KRATOS_CATCH("");
}

template<class TPrimalBaseType>
std::string RansKOmegaNewtonianLaw<TPrimalBaseType>::Info() const
{
    return "RansKOmega" + TPrimalBaseType::Info();
}

template<class TPrimalBaseType>
double RansKOmegaNewtonianLaw<TPrimalBaseType>::GetEffectiveViscosity(ConstitutiveLaw::Parameters& rParameters) const
{
    const Properties& r_prop = rParameters.GetMaterialProperties();
    const double mu = r_prop[DYNAMIC_VISCOSITY];
    const double density = r_prop[DENSITY];

    return mu + density * this->CalculateTurbulentViscosity(rParameters);
}

template<class TPrimalBaseType>
double RansKOmegaNewtonianLaw<TPrimalBaseType>::CalculateTurbulentViscosity(
    ConstitutiveLaw::Parameters& rParameters) const
{
    double tke, omega;
    FluidCalculationUtilities::EvaluateInPoint(
        rParameters.GetElementGeometry(), rParameters.GetShapeFunctionsValues(),
        std::tie(tke, TURBULENT_KINETIC_ENERGY),
        std::tie(omega, TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE));

    double nu_t = 1e-12;
    if (tke > 0.0 && omega > 0.0) {
        nu_t = std::max(tke / omega, 1e-12);
    }

    return nu_t;
}

template<class TPrimalBaseType>
double RansKOmegaNewtonianLaw<TPrimalBaseType>::GetDynamicViscosity(
    ConstitutiveLaw::Parameters& rParameters) const
{
    const Properties& r_prop = rParameters.GetMaterialProperties();
    return r_prop[DYNAMIC_VISCOSITY];
}

template<class TPrimalBaseType>
void RansKOmegaNewtonianLaw<TPrimalBaseType>::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseType)
}

template<class TPrimalBaseType>
void RansKOmegaNewtonianLaw<TPrimalBaseType>::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseType)
}

// template instantiations
template class RansKOmegaNewtonianLaw<Newtonian2DLaw>;
template class RansKOmegaNewtonianLaw<Newtonian3DLaw>;

} // Namespace Kratos
