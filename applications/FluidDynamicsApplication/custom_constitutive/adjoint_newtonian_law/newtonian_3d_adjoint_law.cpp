//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   EffectiveViscositylti-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

// Application includes

// Include base h
#include "newtonian_3d_adjoint_law.h"

namespace Kratos {

// Life cycle /////////////////////////////////////////////////////////////////

Newtonian3DAdjointLaw::Newtonian3DAdjointLaw(ConstitutiveLaw& rConstitutiveLaw):
    BaseType(rConstitutiveLaw) {}

Newtonian3DAdjointLaw::Newtonian3DAdjointLaw(const Newtonian3DAdjointLaw& rOther):
    BaseType(rOther) {}

Newtonian3DAdjointLaw::~Newtonian3DAdjointLaw() {}

// Public operations //////////////////////////////////////////////////////////

FluidAdjointConstitutiveLaw::Pointer Newtonian3DAdjointLaw::Clone() const
{
    return Kratos::make_shared<Newtonian3DAdjointLaw>(*this);
}

void Newtonian3DAdjointLaw::CalculateMaterialResponseCauchyDerivative(
    ConstitutiveLaw::Parameters& rValuesDerivative,
    ConstitutiveLaw::Parameters& rValues,
    const IndexType NodeIndex,
    const Variable<double>& rDerivativeVariable,
    const double EffectiveViscosity,
    const double EffectiveViscosityDerivative)
{
    const Flags& options = rValues.GetOptions();
    const Vector& r_strain_rate = rValues.GetStrainVector();
    const Vector& r_strain_rate_derivative = rValuesDerivative.GetStrainVector();

    const double trace = r_strain_rate[0] + r_strain_rate[1] + r_strain_rate[2];
    const double trace_derivative = r_strain_rate_derivative[0] + r_strain_rate_derivative[1] + r_strain_rate_derivative[2];
    const double volumetric_part = trace / 3.0; // Note: this should be small for an incompressible fluid (it is basically the incompressibility error)
    const double volumetric_part_derivative = trace_derivative / 3.0;

    //computation of stress
    Vector& r_viscous_stress = rValuesDerivative.GetStressVector();
    r_viscous_stress[0] = 2.0 * EffectiveViscosityDerivative * (r_strain_rate[0] - volumetric_part);
    r_viscous_stress[0] += 2.0 * EffectiveViscosity * (r_strain_rate_derivative[0] - volumetric_part_derivative);

    r_viscous_stress[1] = 2.0 * EffectiveViscosityDerivative * (r_strain_rate[1] - volumetric_part);
    r_viscous_stress[1] += 2.0 * EffectiveViscosity * (r_strain_rate_derivative[1] - volumetric_part_derivative);

    r_viscous_stress[2] = 2.0 * EffectiveViscosityDerivative * (r_strain_rate[2] - volumetric_part);
    r_viscous_stress[2] += 2.0 * EffectiveViscosity * (r_strain_rate_derivative[2] - volumetric_part_derivative);

    r_viscous_stress[3] = EffectiveViscosityDerivative * r_strain_rate[3];
    r_viscous_stress[3] += EffectiveViscosity * r_strain_rate_derivative[3];

    r_viscous_stress[4] = EffectiveViscosityDerivative * r_strain_rate[4];
    r_viscous_stress[4] += EffectiveViscosity * r_strain_rate_derivative[4];

    r_viscous_stress[5] = EffectiveViscosityDerivative * r_strain_rate[5];
    r_viscous_stress[5] += EffectiveViscosity * r_strain_rate_derivative[5];

    if(options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
        BaseType::NewtonianConstitutiveMatrixDerivatives3D(EffectiveViscosityDerivative, rValuesDerivative.GetConstitutiveMatrix());
    }

}

int Newtonian3DAdjointLaw::Check(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const ProcessInfo& rCurrentProcessInfo)
{
    return this->GetPrimalConstitutiveLaw().Check(rMaterialProperties, rElementGeometry, rCurrentProcessInfo);
}

double Newtonian3DAdjointLaw::CalculateEffectiveViscosityDerivative(
    ConstitutiveLaw::Parameters& rValuesDerivative,
    ConstitutiveLaw::Parameters& rValues,
    const IndexType NodeIndex,
    const Variable<double>& rDerivativeVariable)
{
    return 0.0;
}

// Info ///////////////////////////////////////////////////////////////////////

std::string Newtonian3DAdjointLaw::Info() const {
    return "Newtonian3DAdjointLaw";
}

void Newtonian3DAdjointLaw::PrintInfo(std::ostream& rOStream) const {
    rOStream << this->Info();
}

void Newtonian3DAdjointLaw::PrintData(std::ostream& rOStream) const {
    rOStream << this->Info();
}

}