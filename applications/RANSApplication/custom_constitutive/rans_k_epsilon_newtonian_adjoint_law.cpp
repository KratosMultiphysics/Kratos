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

// Project includes
#include "includes/checks.h"
#include "includes/variables.h"

// Application includes
#include "custom_utilities/fluid_calculation_utilities.h"
#include "rans_application_variables.h"
#include "custom_constitutive/newtonian_law/newtonian_2d_adjoint_law.h"
#include "custom_constitutive/newtonian_law/newtonian_3d_adjoint_law.h"

// Include base h
#include "rans_k_epsilon_newtonian_adjoint_law.h"

namespace Kratos {

// Life cycle /////////////////////////////////////////////////////////////////

template<class TBaseClassType>
RansKEpsilonNewtonianAdjointLaw<TBaseClassType>::RansKEpsilonNewtonianAdjointLaw(ConstitutiveLaw& rConstitutiveLaw):
    BaseType(rConstitutiveLaw) {}

template<class TBaseClassType>
RansKEpsilonNewtonianAdjointLaw<TBaseClassType>::RansKEpsilonNewtonianAdjointLaw(const RansKEpsilonNewtonianAdjointLaw& rOther):
    BaseType(rOther) {}

template<class TBaseClassType>
RansKEpsilonNewtonianAdjointLaw<TBaseClassType>::~RansKEpsilonNewtonianAdjointLaw() {}

// Public operations //////////////////////////////////////////////////////////

template<class TBaseClassType>
int RansKEpsilonNewtonianAdjointLaw<TBaseClassType>::Check(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    KRATOS_ERROR_IF(rMaterialProperties[DENSITY] <= 0.0)
        << "Incorrect or missing DENSITY provided in material properties "
           "for RansNewtonian2DLaw: "
        << rMaterialProperties[DENSITY] << std::endl;

    KRATOS_ERROR_IF_NOT(rCurrentProcessInfo.Has(TURBULENCE_RANS_C_MU))
        << "TURBULENCE_RANS_C_MU is not defined in the process info.\n";

    for (IndexType i = 0; i < rElementGeometry.PointsNumber(); ++i) {
        const auto& r_node = rElementGeometry[i];

        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TURBULENT_KINETIC_ENERGY, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TURBULENT_ENERGY_DISSIPATION_RATE, r_node);
    }

    return this->GetPrimalConstitutiveLaw().Check(rMaterialProperties, rElementGeometry, rCurrentProcessInfo);

    KRATOS_CATCH("");
}

template<class TBaseClassType>
double RansKEpsilonNewtonianAdjointLaw<TBaseClassType>::CalculateEffectiveViscosityDerivative(
    ConstitutiveLaw::Parameters& rValuesDerivative,
    ConstitutiveLaw::Parameters& rValues,
    const IndexType NodeIndex,
    const Variable<double>& rDerivativeVariable)
{
    const Properties& r_prop = rValues.GetMaterialProperties();
    const double density = r_prop[DENSITY];

    double nu_t_derivative = 0.0;

    // computing gauss point nu_t derivative
    if (rDerivativeVariable == TURBULENT_KINETIC_ENERGY) {
        nu_t_derivative = rValues.GetElementGeometry().GetValue(TURBULENT_VISCOSITY_DERIVATIVES)[0];
    } else if (rDerivativeVariable == TURBULENT_ENERGY_DISSIPATION_RATE) {
        nu_t_derivative = rValues.GetElementGeometry().GetValue(TURBULENT_VISCOSITY_DERIVATIVES)[1];
    }

    return density * nu_t_derivative;
}

// Info ///////////////////////////////////////////////////////////////////////

template<class TBaseClassType>
std::string RansKEpsilonNewtonianAdjointLaw<TBaseClassType>::Info() const {
    return "RansKEpsilonNewtonianAdjointLaw";
}

template<class TBaseClassType>
void RansKEpsilonNewtonianAdjointLaw<TBaseClassType>::PrintInfo(std::ostream& rOStream) const {
    rOStream << this->Info();
}

template<class TBaseClassType>
void RansKEpsilonNewtonianAdjointLaw<TBaseClassType>::PrintData(std::ostream& rOStream) const {
    rOStream << this->Info();
}

// template instantiations

template class RansKEpsilonNewtonianAdjointLaw<Newtonian2DAdjointLaw>;
template class RansKEpsilonNewtonianAdjointLaw<Newtonian3DAdjointLaw>;

}