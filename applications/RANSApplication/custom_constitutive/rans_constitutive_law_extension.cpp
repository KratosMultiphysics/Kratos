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
#include <sstream>

// External includes

// Project includes
#include "includes/cfd_variables.h"
#include "includes/checks.h"

// Application includes
#include "custom_utilities/fluid_calculation_utilities.h"

#include "custom_constitutive/newtonian_2d_law.h"
#include "custom_constitutive/newtonian_3d_law.h"

// Include base h
#include "rans_constitutive_law_extension.h"

namespace Kratos
{
///@name Classes
///@{

template <class TFluidConstitutiveLawType>
RansConstitutiveLawExtension<TFluidConstitutiveLawType>::RansConstitutiveLawExtension()
    : BaseType()
{
}

template <class TFluidConstitutiveLawType>
RansConstitutiveLawExtension<TFluidConstitutiveLawType>::RansConstitutiveLawExtension(
    const RansConstitutiveLawExtension& rOther)
    : BaseType(rOther)
{
}

template <class TFluidConstitutiveLawType>
ConstitutiveLaw::Pointer RansConstitutiveLawExtension<TFluidConstitutiveLawType>::Clone() const
{
    return Kratos::make_shared<RansConstitutiveLawExtension<TFluidConstitutiveLawType>>(*this);
}

template <class TFluidConstitutiveLawType>
RansConstitutiveLawExtension<TFluidConstitutiveLawType>::~RansConstitutiveLawExtension()
{
}

template <class TFluidConstitutiveLawType>
int RansConstitutiveLawExtension<TFluidConstitutiveLawType>::Check(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    int value = BaseType::Check(rMaterialProperties, rElementGeometry, rCurrentProcessInfo);

    for (IndexType i = 0; i < rElementGeometry.PointsNumber(); ++i) {
        const auto& r_node = rElementGeometry[i];

        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TURBULENT_VISCOSITY, r_node);
    }

    return value;

    KRATOS_CATCH("");
}

template <class TFluidConstitutiveLawType>
std::string RansConstitutiveLawExtension<TFluidConstitutiveLawType>::Info() const
{
    std::stringstream msg;
    msg << "Rans" << BaseType::Info();
    return msg.str();
}

template <class TFluidConstitutiveLawType>
double RansConstitutiveLawExtension<TFluidConstitutiveLawType>::GetEffectiveViscosity(
    ConstitutiveLaw::Parameters& rParameters) const
{
    const double mu = BaseType::GetEffectiveViscosity(rParameters);
    const double density = rParameters.GetMaterialProperties().GetValue(DENSITY);

    double turbulent_nu;
    FluidCalculationUtilities::EvaluateInPoint(
        rParameters.GetElementGeometry(), rParameters.GetShapeFunctionsValues(),
        std::tie(turbulent_nu, TURBULENT_VISCOSITY));

    return mu + density * turbulent_nu;
}

template <class TFluidConstitutiveLawType>
void RansConstitutiveLawExtension<TFluidConstitutiveLawType>::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseType);
}

template <class TFluidConstitutiveLawType>
void RansConstitutiveLawExtension<TFluidConstitutiveLawType>::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseType);
}

// template instantiations
template class RansConstitutiveLawExtension<Newtonian2DLaw>;
template class RansConstitutiveLawExtension<Newtonian3DLaw>;

} // namespace Kratos.
