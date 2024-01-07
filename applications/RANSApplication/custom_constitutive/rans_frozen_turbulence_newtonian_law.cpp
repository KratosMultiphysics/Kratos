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

#include "rans_k_epsilon_newtonian_law.h"
#include "rans_k_omega_newtonian_law.h"
#include "rans_k_omega_sst_newtonian_law.h"

// Include base h
#include "custom_constitutive/rans_frozen_turbulence_newtonian_law.h"

namespace Kratos
{
//******************************CONSTRUCTOR*******************************************
//************************************************************************************

template<class TUnfrozenTurbulenceConsititutiveLaw>
RansFrozenTurbulenceNewtonianLaw<TUnfrozenTurbulenceConsititutiveLaw>::RansFrozenTurbulenceNewtonianLaw() : BaseType()
{
}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

template<class TUnfrozenTurbulenceConsititutiveLaw>
RansFrozenTurbulenceNewtonianLaw<TUnfrozenTurbulenceConsititutiveLaw>::RansFrozenTurbulenceNewtonianLaw(const RansFrozenTurbulenceNewtonianLaw& rOther)
    : BaseType(rOther)
{
}

//********************************CLONE***********************************************
//************************************************************************************

template<class TUnfrozenTurbulenceConsititutiveLaw>
ConstitutiveLaw::Pointer RansFrozenTurbulenceNewtonianLaw<TUnfrozenTurbulenceConsititutiveLaw>::Clone() const
{
    return Kratos::make_shared<RansFrozenTurbulenceNewtonianLaw>(*this);
}

template<class TUnfrozenTurbulenceConsititutiveLaw>
void RansFrozenTurbulenceNewtonianLaw<TUnfrozenTurbulenceConsititutiveLaw>::CalculateDerivative(
    ConstitutiveLaw::Parameters& rParameters,
    const Variable<double>& rFunctionVariable,
    const Variable<double>& rDerivativeVariable,
    double& rOutput)
{
    KRATOS_TRY

    if (rFunctionVariable == EFFECTIVE_VISCOSITY) {
        rOutput = 0.0;
    }

    KRATOS_CATCH("");
}

template<class TUnfrozenTurbulenceConsititutiveLaw>
std::string RansFrozenTurbulenceNewtonianLaw<TUnfrozenTurbulenceConsititutiveLaw>::Info() const
{
    return "RansFrozenTurbulence" + TUnfrozenTurbulenceConsititutiveLaw::Info();
}

template<class TUnfrozenTurbulenceConsititutiveLaw>
void RansFrozenTurbulenceNewtonianLaw<TUnfrozenTurbulenceConsititutiveLaw>::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseType)
}

template<class TUnfrozenTurbulenceConsititutiveLaw>
void RansFrozenTurbulenceNewtonianLaw<TUnfrozenTurbulenceConsititutiveLaw>::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseType)
}

// template instantiations
template class RansFrozenTurbulenceNewtonianLaw<RansKEpsilonNewtonianLaw<Newtonian2DLaw>>;
template class RansFrozenTurbulenceNewtonianLaw<RansKEpsilonNewtonianLaw<Newtonian3DLaw>>;

template class RansFrozenTurbulenceNewtonianLaw<RansKOmegaNewtonianLaw<Newtonian2DLaw>>;
template class RansFrozenTurbulenceNewtonianLaw<RansKOmegaNewtonianLaw<Newtonian3DLaw>>;

template class RansFrozenTurbulenceNewtonianLaw<RansKOmegaSSTNewtonianLaw<2, Newtonian2DLaw>>;
template class RansFrozenTurbulenceNewtonianLaw<RansKOmegaSSTNewtonianLaw<3, Newtonian3DLaw>>;

} // Namespace Kratos
