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

// Application includes

// Include base h
#include "newtonian_2d_adjoint_law.h"

namespace Kratos {

// Life cycle /////////////////////////////////////////////////////////////////

Newtonian2DAdjointLaw::Newtonian2DAdjointLaw(ConstitutiveLaw& rConstitutiveLaw):
    BaseType(rConstitutiveLaw) {}

Newtonian2DAdjointLaw::Newtonian2DAdjointLaw(const Newtonian2DAdjointLaw& rOther):
    BaseType(rOther) {}

Newtonian2DAdjointLaw::~Newtonian2DAdjointLaw() {}

// Public operations //////////////////////////////////////////////////////////

void Newtonian2DAdjointLaw::CalculateMaterialResponseCauchyDerivative(
    ConstitutiveLaw::Parameters& rValuesDerivative,
    ConstitutiveLaw::Parameters& rValues,
    const IndexType NodeIndex,
    const Variable<double>& rDerivativeVariable,
    const double EffectiveViscosity,
    const double EffectiveViscosityDerivative)
{
    const Flags& options = rValues.GetOptions();
    // const Vector& r_strain_rate = rValues.GetStrainVector();
    // const Vector& r_strain_rate_derivative = rValuesDerivative.GetStrainVector();

    // const double trace = r_strain_rate[0] + r_strain_rate[1];
    // const double trace_derivative = r_strain_rate_derivative[0] + r_strain_rate_derivative[1];
    // const double volumetric_part = trace / 3.0; // Note: this should be small for an incompressible fluid (it is basically the incompressibility error)
    // const double volumetric_part_derivative = trace_derivative / 3.0;

    // //computation of stress
    // Vector& r_viscous_stress_derivative = rValuesDerivative.GetStressVector();

    // r_viscous_stress_derivative[0] = 2.0 * EffectiveViscosityDerivative * (r_strain_rate[0] - volumetric_part);
    // r_viscous_stress_derivative[0] += 2.0 * EffectiveViscosity * (r_strain_rate_derivative[0] - volumetric_part_derivative);

    // r_viscous_stress_derivative[1] = 2.0 * EffectiveViscosityDerivative * (r_strain_rate[1] - volumetric_part);
    // r_viscous_stress_derivative[1] += 2.0 * EffectiveViscosity * (r_strain_rate_derivative[1] - volumetric_part_derivative);

    // r_viscous_stress_derivative[2] = EffectiveViscosityDerivative * r_strain_rate[2];
    // r_viscous_stress_derivative[2] += EffectiveViscosity * r_strain_rate_derivative[2];

    // if(options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
    //     BaseType::NewtonianConstitutiveMatrixDerivatives2D(EffectiveViscosityDerivative, rValuesDerivative.GetConstitutiveMatrix());
    // }

}

int Newtonian2DAdjointLaw::Check(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const ProcessInfo& rCurrentProcessInfo)
{
    return this->GetPrimalConstitutiveLaw().Check(rMaterialProperties, rElementGeometry, rCurrentProcessInfo);
}

double Newtonian2DAdjointLaw::CalculateEffectiveViscosityDerivative(
    ConstitutiveLaw::Parameters& rValuesDerivative,
    ConstitutiveLaw::Parameters& rValues,
    const IndexType NodeIndex,
    const Variable<double>& rDerivativeVariable)
{
    return 0.0;
}

// Info ///////////////////////////////////////////////////////////////////////

std::string Newtonian2DAdjointLaw::Info() const {
    return "Newtonian2DAdjointLaw";
}

void Newtonian2DAdjointLaw::PrintInfo(std::ostream& rOStream) const {
    rOStream << this->Info();
}

void Newtonian2DAdjointLaw::PrintData(std::ostream& rOStream) const {
    rOStream << this->Info();
}

}