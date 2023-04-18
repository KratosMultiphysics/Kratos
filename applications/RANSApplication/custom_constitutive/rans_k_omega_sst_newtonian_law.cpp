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
#include "custom_elements/data_containers/k_omega_sst/element_data_utilities.h"
#include "custom_utilities/fluid_calculation_utilities.h"
#include "rans_application_variables.h"

// Include base h
#include "custom_constitutive/rans_k_omega_sst_newtonian_law.h"

namespace Kratos
{
//******************************CONSTRUCTOR*******************************************
//************************************************************************************

template<unsigned int TDim, class TPrimalBaseType>
RansKOmegaSSTNewtonianLaw<TDim, TPrimalBaseType>::RansKOmegaSSTNewtonianLaw() : BaseType()
{
}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

template<unsigned int TDim, class TPrimalBaseType>
RansKOmegaSSTNewtonianLaw<TDim, TPrimalBaseType>::RansKOmegaSSTNewtonianLaw(const RansKOmegaSSTNewtonianLaw& rOther)
    : BaseType(rOther)
{
}

//********************************CLONE***********************************************
//************************************************************************************

template<unsigned int TDim, class TPrimalBaseType>
ConstitutiveLaw::Pointer RansKOmegaSSTNewtonianLaw<TDim, TPrimalBaseType>::Clone() const
{
    return Kratos::make_shared<RansKOmegaSSTNewtonianLaw>(*this);
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

template<unsigned int TDim, class TPrimalBaseType>
RansKOmegaSSTNewtonianLaw<TDim, TPrimalBaseType>::~RansKOmegaSSTNewtonianLaw()
{
}

template<unsigned int TDim, class TPrimalBaseType>
int RansKOmegaSSTNewtonianLaw<TDim, TPrimalBaseType>::Check(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    // Check viscosity value
    KRATOS_ERROR_IF(rMaterialProperties[DYNAMIC_VISCOSITY] <= 0.0)
        << "Incorrect or missing DYNAMIC_VISCOSITY provided in material properties "
           "for RansKOmegaSSTNewtonianLaw: "
        << rMaterialProperties[DYNAMIC_VISCOSITY] << std::endl;

    KRATOS_ERROR_IF(rMaterialProperties[DENSITY] <= 0.0)
        << "Incorrect or missing DENSITY provided in material properties "
           "for RansKOmegaSSTNewtonianLaw: "
        << rMaterialProperties[DENSITY] << std::endl;

    KRATOS_ERROR_IF_NOT(rCurrentProcessInfo.Has(TURBULENCE_RANS_C_MU))
        << "TURBULENCE_RANS_C_MU is not found in process info.\n";

    KRATOS_ERROR_IF_NOT(rCurrentProcessInfo.Has(TURBULENCE_RANS_A1))
        << "TURBULENCE_RANS_A1 is not found in process info.\n";

    return 0;

    KRATOS_CATCH("");
}

template<unsigned int TDim, class TPrimalBaseType>
double& RansKOmegaSSTNewtonianLaw<TDim, TPrimalBaseType>::CalculateValue(
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

template<unsigned int TDim, class TPrimalBaseType>
std::string RansKOmegaSSTNewtonianLaw<TDim, TPrimalBaseType>::Info() const
{
    return "RansKOmegaSST" + TPrimalBaseType::Info();
}

template<unsigned int TDim, class TPrimalBaseType>
double RansKOmegaSSTNewtonianLaw<TDim, TPrimalBaseType>::GetEffectiveViscosity(ConstitutiveLaw::Parameters& rParameters) const
{
    const Properties& r_prop = rParameters.GetMaterialProperties();
    const double mu = r_prop[DYNAMIC_VISCOSITY];
    const double density = r_prop[DENSITY];

    return mu + density * this->CalculateTurbulentViscosity(rParameters);
}

template<unsigned int TDim, class TPrimalBaseType>
double RansKOmegaSSTNewtonianLaw<TDim, TPrimalBaseType>::CalculateTurbulentViscosity(
    ConstitutiveLaw::Parameters& rParameters) const
{
    const auto& r_process_info = rParameters.GetProcessInfo();
    const double beta_star = r_process_info[TURBULENCE_RANS_C_MU];
    const double a1 = r_process_info[TURBULENCE_RANS_A1];

    const auto& r_properties = rParameters.GetMaterialProperties();
    const double rho = r_properties[DENSITY];
    const double nu = r_properties[DYNAMIC_VISCOSITY] / rho;

    double tke, omega, y;
    FluidCalculationUtilities::EvaluateInPoint(
        rParameters.GetElementGeometry(), rParameters.GetShapeFunctionsValues(),
        std::tie(tke, TURBULENT_KINETIC_ENERGY),
        std::tie(omega, TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE),
        std::tie(y, DISTANCE)
    );

    BoundedMatrix<double, TDim, TDim> velocity_gradient;
    FluidCalculationUtilities::EvaluateGradientInPoint(
        rParameters.GetElementGeometry(), rParameters.GetShapeFunctionsDerivatives(),
        std::tie(velocity_gradient, VELOCITY));

    const double f_2 = KOmegaSSTElementData::CalculateF2(tke, omega, nu, y, beta_star);

    const BoundedMatrix<double, TDim, TDim> symmetric_velocity_gradient =
        (velocity_gradient + trans(velocity_gradient)) * 0.5;

    const double t = norm_frobenius(symmetric_velocity_gradient) * 1.414;

    return std::max(KOmegaSSTElementData::CalculateTurbulentKinematicViscosity(tke, omega, t, f_2, a1), 1e-12);
}

template<unsigned int TDim, class TPrimalBaseType>
double RansKOmegaSSTNewtonianLaw<TDim, TPrimalBaseType>::GetDynamicViscosity(
    ConstitutiveLaw::Parameters& rParameters) const
{
    const Properties& r_prop = rParameters.GetMaterialProperties();
    return r_prop[DYNAMIC_VISCOSITY];
}

template<unsigned int TDim, class TPrimalBaseType>
void RansKOmegaSSTNewtonianLaw<TDim, TPrimalBaseType>::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseType)
}

template<unsigned int TDim, class TPrimalBaseType>
void RansKOmegaSSTNewtonianLaw<TDim, TPrimalBaseType>::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseType)
}

// template instantiations
template class RansKOmegaSSTNewtonianLaw<2, Newtonian2DLaw>;
template class RansKOmegaSSTNewtonianLaw<3, Newtonian3DLaw>;

} // Namespace Kratos
