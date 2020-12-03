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

// Application includes
#include "custom_constitutive/fluid_constitutive_law.h"

// Include base h
#include "fluid_adjoint_constitutive_law.h"

namespace Kratos {

// Life cycle /////////////////////////////////////////////////////////////////

FluidAdjointConstitutiveLaw::FluidAdjointConstitutiveLaw(ConstitutiveLaw& rConstitutiveLaw):
    mrPrimalConstitutiveLaw(rConstitutiveLaw) {}

FluidAdjointConstitutiveLaw::FluidAdjointConstitutiveLaw(const FluidAdjointConstitutiveLaw& rOther):
    mrPrimalConstitutiveLaw(rOther.mrPrimalConstitutiveLaw) {}

FluidAdjointConstitutiveLaw::~FluidAdjointConstitutiveLaw() {}

// Public operations //////////////////////////////////////////////////////////

void FluidAdjointConstitutiveLaw::CalculateMaterialResponseCauchyDerivative(
    ConstitutiveLaw::Parameters& rValuesDerivative,
    ConstitutiveLaw::Parameters& rValues,
    const IndexType NodeIndex,
    const Variable<double>& rDerivativeVariable,
    const double EffectiveViscosity,
    const double EffectiveViscosityDerivative)
{
    KRATOS_ERROR << "Calling base "
                    "FluidAdjointConstitutiveLaw::CalculateMaterialResponseCauchyDerivative "
                    "method. This class should not be instantiated. Please "
                    "check your constitutive law."
                 << std::endl;
}

int FluidAdjointConstitutiveLaw::Check(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_ERROR << "Calling base "
                    "FluidAdjointConstitutiveLaw::Check "
                    "method. This class should not be instantiated. Please "
                    "check your constitutive law."
                 << std::endl;
    return 999;
}

double FluidAdjointConstitutiveLaw::CalculateEffectiveViscosityDerivative(
    ConstitutiveLaw::Parameters& rValuesDerivative,
    ConstitutiveLaw::Parameters& rValues,
    const IndexType NodeIndex,
    const Variable<double>& rDerivativeVariable)
{
    KRATOS_ERROR << "Calling base "
                    "FluidAdjointConstitutiveLaw::GetEffectiveViscosityDerivative "
                    "method. This class should not be instantiated. Please "
                    "check your constitutive law."
                 << std::endl;
    return 0.0;
}

// Inquiry ////////////////////////////////////////////////////////////////////

ConstitutiveLaw& FluidAdjointConstitutiveLaw::GetPrimalConstitutiveLaw()
{
    return mrPrimalConstitutiveLaw;
}

// Info ///////////////////////////////////////////////////////////////////////

std::string FluidAdjointConstitutiveLaw::Info() const {
    return "FluidAdjointConstitutiveLaw";
}

void FluidAdjointConstitutiveLaw::PrintInfo(std::ostream& rOStream) const {
    rOStream << this->Info();
}

void FluidAdjointConstitutiveLaw::PrintData(std::ostream& rOStream) const {
    rOStream << this->Info();
}

// Protected operations //////////////////////////////////////////////////////

void FluidAdjointConstitutiveLaw::NewtonianConstitutiveMatrixDerivatives2D(
    const double EffectiveViscosityDerivative,
    Matrix& rCDerivative) const
{
    FluidConstitutiveLaw::NewtonianConstitutiveMatrix2D(EffectiveViscosityDerivative, rCDerivative);
}

void FluidAdjointConstitutiveLaw::NewtonianConstitutiveMatrixDerivatives3D(
    const double EffectiveViscosityDerivative,
    Matrix& rCDerivative) const
{
    FluidConstitutiveLaw::NewtonianConstitutiveMatrix3D(EffectiveViscosityDerivative, rCDerivative);
}

}