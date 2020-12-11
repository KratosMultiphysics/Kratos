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
#include "custom_utilities/fluid_calculation_utilities.h"

// Include base h
#include "custom_constitutive/rans_newtonian_2d_law.h"

namespace Kratos
{
//******************************CONSTRUCTOR*******************************************
//************************************************************************************

RansNewtonian2DLaw::RansNewtonian2DLaw() : BaseType()
{
}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

RansNewtonian2DLaw::RansNewtonian2DLaw(const RansNewtonian2DLaw& rOther)
    : BaseType(rOther)
{
}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer RansNewtonian2DLaw::Clone() const
{
    return Kratos::make_shared<RansNewtonian2DLaw>(*this);
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

RansNewtonian2DLaw::~RansNewtonian2DLaw()
{
}

ConstitutiveLaw::SizeType RansNewtonian2DLaw::WorkingSpaceDimension()
{
    return 2;
}

ConstitutiveLaw::SizeType RansNewtonian2DLaw::GetStrainSize()
{
    return 3;
}

void RansNewtonian2DLaw::CalculateMaterialResponseCauchy(Parameters& rValues)
{
    const Flags& options = rValues.GetOptions();
    const Vector& r_strain_rate = rValues.GetStrainVector();
    Vector& r_viscous_stress = rValues.GetStressVector();

    const double mu = this->GetEffectiveViscosity(rValues);

    const double trace = r_strain_rate[0] + r_strain_rate[1];
    const double volumetric_part =
        trace / 3.0; // Note: this should be small for an incompressible fluid (it is basically the incompressibility error)

    // computation of stress
    r_viscous_stress[0] = 2.0 * mu * (r_strain_rate[0] - volumetric_part);
    r_viscous_stress[1] = 2.0 * mu * (r_strain_rate[1] - volumetric_part);
    r_viscous_stress[2] = mu * r_strain_rate[2];

    if (options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
        this->NewtonianConstitutiveMatrix2D(mu, rValues.GetConstitutiveMatrix());
    }
}

int RansNewtonian2DLaw::Check(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Check viscosity value
    KRATOS_ERROR_IF(rMaterialProperties[DYNAMIC_VISCOSITY] <= 0.0)
        << "Incorrect or missing DYNAMIC_VISCOSITY provided in process info "
           "for RansNewtonian2DLaw: "
        << rMaterialProperties[DYNAMIC_VISCOSITY] << std::endl;

    KRATOS_ERROR_IF(rMaterialProperties[DENSITY] <= 0.0)
        << "Incorrect or missing DENSITY provided in process info "
           "for RansNewtonian2DLaw: "
        << rMaterialProperties[DENSITY] << std::endl;

    for (IndexType i = 0; i < rElementGeometry.PointsNumber(); ++i) {
        const auto& r_node = rElementGeometry[i];

        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TURBULENT_VISCOSITY, r_node);
    }

    return 0;

    KRATOS_CATCH("");
}

std::string RansNewtonian2DLaw::Info() const
{
    return "RansNewtonian2DLaw";
}

double RansNewtonian2DLaw::GetEffectiveViscosity(ConstitutiveLaw::Parameters& rParameters) const
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

void RansNewtonian2DLaw::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, FluidConstitutiveLaw)
}

void RansNewtonian2DLaw::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, FluidConstitutiveLaw)
}

} // Namespace Kratos
