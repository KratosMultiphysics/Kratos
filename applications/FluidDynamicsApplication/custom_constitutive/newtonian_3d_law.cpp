//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//

// System includes
#include <iostream>

// External includes

// Project includes
#include "includes/cfd_variables.h"
#include "includes/checks.h"

// Application includes
#include "fluid_dynamics_application_variables.h"

// Include base h
#include "custom_constitutive/newtonian_3d_law.h"

namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

Newtonian3DLaw::Newtonian3DLaw()
    : FluidConstitutiveLaw()
{}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

Newtonian3DLaw::Newtonian3DLaw(const Newtonian3DLaw& rOther)
    : FluidConstitutiveLaw(rOther)
{}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer Newtonian3DLaw::Clone() const {
    return Kratos::make_shared<Newtonian3DLaw>(*this);
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

Newtonian3DLaw::~Newtonian3DLaw() {}

ConstitutiveLaw::SizeType Newtonian3DLaw::WorkingSpaceDimension() {
    return 3;
}

ConstitutiveLaw::SizeType Newtonian3DLaw::GetStrainSize() const {
    return 6;
}

void  Newtonian3DLaw::CalculateMaterialResponseCauchy (Parameters& rValues)
{
    const Flags& options = rValues.GetOptions();
    const Vector& r_strain_rate = rValues.GetStrainVector();
    Vector& r_viscous_stress = rValues.GetStressVector();

    const double mu = this->GetEffectiveViscosity(rValues);

    const double trace = r_strain_rate[0] + r_strain_rate[1] + r_strain_rate[2];
    const double volumetric_part = trace/3.0; // Note: this should be small for an incompressible fluid (it is basically the incompressibility error)

    //computation of stress
    r_viscous_stress[0] = 2.0*mu*(r_strain_rate[0] - volumetric_part);
    r_viscous_stress[1] = 2.0*mu*(r_strain_rate[1] - volumetric_part);
    r_viscous_stress[2] = 2.0*mu*(r_strain_rate[2] - volumetric_part);
    r_viscous_stress[3] = mu*r_strain_rate[3];
    r_viscous_stress[4] = mu*r_strain_rate[4];
    r_viscous_stress[5] = mu*r_strain_rate[5];

    if( options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ) )
    {
        this->NewtonianConstitutiveMatrix3D(mu,rValues.GetConstitutiveMatrix());
    }
}

void Newtonian3DLaw::CalculateDerivative(
    Parameters& rParameterValues,
    const Variable<double>& rFunctionVariable,
    const Variable<double>& rDerivativeVariable,
    double& rOutput)
{
    if (rFunctionVariable == EFFECTIVE_VISCOSITY) {
        // since EFFECTIVE_VISCOSITY is a constant for
        // fluid domain, this derivative remains zero.
        rOutput = 0.0;
    } else {
        FluidConstitutiveLaw::CalculateDerivative(
            rParameterValues, rFunctionVariable, rDerivativeVariable, rOutput);
    }
}

void Newtonian3DLaw::CalculateDerivative(
    Parameters& rParameterValues,
    const Variable<Vector>& rFunctionVariable,
    const Variable<double>& rDerivativeVariable,
    Vector& rOutput)
{
    KRATOS_TRY

    if (rFunctionVariable == CAUCHY_STRESS_VECTOR) {
        // computes derivatives of CAUCHY_STRESS_VECTOR
        // resize and clear the output variable
        if (rOutput.size() != 6) {
            rOutput.resize(6);
        }
        rOutput.clear();

        if (rDerivativeVariable.IsComponent() && rDerivativeVariable.GetSourceVariable() == STRAIN_RATE_3D) {
            // compute derivatives w.r.t. STRAIN_RATE
            const double mu = this->GetEffectiveViscosity(rParameterValues);
            const double volumetric_part_derivative = (rDerivativeVariable.GetComponentIndex() < 3) ? 1.0 / 3.0 : 0.0;

            rOutput[0] = 2.0 * mu * ((rDerivativeVariable.GetComponentIndex() == 0) - volumetric_part_derivative);
            rOutput[1] = 2.0 * mu * ((rDerivativeVariable.GetComponentIndex() == 1) - volumetric_part_derivative);
            rOutput[2] = 2.0 * mu * ((rDerivativeVariable.GetComponentIndex() == 2) - volumetric_part_derivative);
            rOutput[3] = mu * (rDerivativeVariable.GetComponentIndex() == 3);
            rOutput[4] = mu * (rDerivativeVariable.GetComponentIndex() == 4);
            rOutput[5] = mu * (rDerivativeVariable.GetComponentIndex() == 5);

        } else if (rDerivativeVariable == EFFECTIVE_VISCOSITY) {
            // compute derivatives w.r.t. effective viscosity
            const Vector& r_strain_rate = rParameterValues.GetStrainVector();
            const double trace = r_strain_rate[0] + r_strain_rate[1] + r_strain_rate[2];
            const double volumetric_part = trace / 3.0;

            rOutput[0] = 2.0 * (r_strain_rate[0] - volumetric_part);
            rOutput[1] = 2.0 * (r_strain_rate[1] - volumetric_part);
            rOutput[2] = 2.0 * (r_strain_rate[2] - volumetric_part);
            rOutput[3] = r_strain_rate[3];
            rOutput[4] = r_strain_rate[4];
            rOutput[5] = r_strain_rate[5];
        } else {
            FluidConstitutiveLaw::CalculateDerivative(
                rParameterValues, rFunctionVariable, rDerivativeVariable, rOutput);
        }
    } else {
        FluidConstitutiveLaw::CalculateDerivative(
            rParameterValues, rFunctionVariable, rDerivativeVariable, rOutput);
    }

    KRATOS_CATCH("");
}

void Newtonian3DLaw::CalculateDerivative(
    Parameters& rParameterValues,
    const Variable<Matrix>& rFunctionVariable,
    const Variable<double>& rDerivativeVariable,
    Matrix& rOutput)
{
    KRATOS_TRY

    if (rFunctionVariable == CONSTITUTIVE_MATRIX) {
        // computes derivatives of CONSTITUTIVE_MATRIX
        // resize and clear the output variable
        if (rOutput.size1() != 6 || rOutput.size2() != 6) {
            rOutput.resize(6, 6, false);
        }
        rOutput.clear();

        if (rDerivativeVariable == EFFECTIVE_VISCOSITY) {
            // compute derivatives w.r.t. effective viscosity
            FluidConstitutiveLaw::NewtonianConstitutiveMatrix3D(1.0, rOutput);
        } else {
            FluidConstitutiveLaw::CalculateDerivative(
                rParameterValues, rFunctionVariable, rDerivativeVariable, rOutput);
        }
    } else {
        FluidConstitutiveLaw::CalculateDerivative(
            rParameterValues, rFunctionVariable, rDerivativeVariable, rOutput);
    }

    KRATOS_CATCH("");
}

int Newtonian3DLaw::Check(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const ProcessInfo& rCurrentProcessInfo) const
{
    // Check viscosity value
    KRATOS_ERROR_IF(rMaterialProperties[DYNAMIC_VISCOSITY] <= 0.0)
        << "Incorrect or missing DYNAMIC_VISCOSITY provided in process info for Newtonian3DLaw: " << rMaterialProperties[DYNAMIC_VISCOSITY] << std::endl;

    return 0;
}

std::string Newtonian3DLaw::Info() const {
    return "Newtonian3DLaw";
}

double Newtonian3DLaw::GetEffectiveViscosity(ConstitutiveLaw::Parameters& rParameters) const
{
    const Properties &r_prop = rParameters.GetMaterialProperties();
    const double effective_viscosity = r_prop[DYNAMIC_VISCOSITY];
    return effective_viscosity;
}

void Newtonian3DLaw::save(Serializer& rSerializer) const {
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, FluidConstitutiveLaw )
}

void Newtonian3DLaw::load(Serializer& rSerializer) {
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, FluidConstitutiveLaw )
}

} // Namespace Kratos
