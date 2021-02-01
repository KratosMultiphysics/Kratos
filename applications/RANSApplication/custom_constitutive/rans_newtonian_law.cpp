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
#include "custom_constitutive/newtonian_law/newtonian_2d_adjoint_law.h"
#include "custom_constitutive/newtonian_law/newtonian_3d_adjoint_law.h"
#include "custom_constitutive/rans_k_epsilon_newtonian_adjoint_law.h"
#include "custom_utilities/fluid_calculation_utilities.h"

// Include base h
#include "custom_constitutive/rans_newtonian_law.h"

namespace Kratos
{
//******************************CONSTRUCTOR*******************************************
//************************************************************************************

template<class TPrimalBaseType, class TAdjointBaseType>
RansNewtonianLaw<TPrimalBaseType, TAdjointBaseType>::RansNewtonianLaw() : BaseType()
{
}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

template<class TPrimalBaseType, class TAdjointBaseType>
RansNewtonianLaw<TPrimalBaseType, TAdjointBaseType>::RansNewtonianLaw(const RansNewtonianLaw& rOther)
    : BaseType(rOther)
{
}

//********************************CLONE***********************************************
//************************************************************************************

template<class TPrimalBaseType, class TAdjointBaseType>
ConstitutiveLaw::Pointer RansNewtonianLaw<TPrimalBaseType, TAdjointBaseType>::Clone() const
{
    return Kratos::make_shared<RansNewtonianLaw>(*this);
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

template<class TPrimalBaseType, class TAdjointBaseType>
RansNewtonianLaw<TPrimalBaseType, TAdjointBaseType>::~RansNewtonianLaw()
{
}

template<class TPrimalBaseType, class TAdjointBaseType>
int RansNewtonianLaw<TPrimalBaseType, TAdjointBaseType>::Check(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Check viscosity value
    KRATOS_ERROR_IF(rMaterialProperties[DYNAMIC_VISCOSITY] <= 0.0)
        << "Incorrect or missing DYNAMIC_VISCOSITY provided in material properties "
           "for RansNewtonianLaw: "
        << rMaterialProperties[DYNAMIC_VISCOSITY] << std::endl;

    KRATOS_ERROR_IF(rMaterialProperties[DENSITY] <= 0.0)
        << "Incorrect or missing DENSITY provided in material properties "
           "for RansNewtonianLaw: "
        << rMaterialProperties[DENSITY] << std::endl;

    for (IndexType i = 0; i < rElementGeometry.PointsNumber(); ++i) {
        const auto& r_node = rElementGeometry[i];

        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TURBULENT_VISCOSITY, r_node);
    }

    return 0;

    KRATOS_CATCH("");
}

template<class TPrimalBaseType, class TAdjointBaseType>
FluidAdjointConstitutiveLaw::Pointer RansNewtonianLaw<TPrimalBaseType, TAdjointBaseType>::GetAdjointConstitutiveLaw()
{
    return Kratos::make_shared<TAdjointBaseType>(*this);
}

template<class TPrimalBaseType, class TAdjointBaseType>
std::string RansNewtonianLaw<TPrimalBaseType, TAdjointBaseType>::Info() const
{
    return "Rans" + TPrimalBaseType::Info();
}

template<class TPrimalBaseType, class TAdjointBaseType>
double RansNewtonianLaw<TPrimalBaseType, TAdjointBaseType>::GetEffectiveViscosity(ConstitutiveLaw::Parameters& rParameters) const
{
    const Properties& r_prop = rParameters.GetMaterialProperties();

    const double mu = r_prop[DYNAMIC_VISCOSITY];
    const double density = r_prop[DENSITY];

    double turbulent_nu;
    FluidCalculationUtilities::EvaluateInPoint(
        rParameters.GetElementGeometry(), rParameters.GetShapeFunctionsValues(),
        std::tie(turbulent_nu, TURBULENT_VISCOSITY));

    return mu + density * turbulent_nu;
}

template<class TPrimalBaseType, class TAdjointBaseType>
void RansNewtonianLaw<TPrimalBaseType, TAdjointBaseType>::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseType)
}

template<class TPrimalBaseType, class TAdjointBaseType>
void RansNewtonianLaw<TPrimalBaseType, TAdjointBaseType>::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseType)
}

// template instantiations
template class RansNewtonianLaw<Newtonian2DLaw, RansKEpsilonNewtonianAdjointLaw<Newtonian2DAdjointLaw>>;
template class RansNewtonianLaw<Newtonian3DLaw, RansKEpsilonNewtonianAdjointLaw<Newtonian3DAdjointLaw>>;

} // Namespace Kratos
