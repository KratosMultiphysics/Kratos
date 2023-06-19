//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela
//

#include "fluid_constitutive_law.h"

namespace Kratos {

// Life cycle /////////////////////////////////////////////////////////////////

FluidConstitutiveLaw::FluidConstitutiveLaw():
    ConstitutiveLaw() {}

FluidConstitutiveLaw::FluidConstitutiveLaw(const FluidConstitutiveLaw& rOther):
    ConstitutiveLaw(rOther) {}

FluidConstitutiveLaw::~FluidConstitutiveLaw() {}

// Public operations //////////////////////////////////////////////////////////

ConstitutiveLaw::Pointer FluidConstitutiveLaw::Clone() const {
    KRATOS_ERROR << "Calling base FluidConstitutiveLaw::Clone method. This "
                    "class should not be instantiated. Please check your "
                    "constitutive law."
                 << std::endl;
    return Kratos::make_shared<FluidConstitutiveLaw>(*this);
}

void FluidConstitutiveLaw::CalculateMaterialResponseCauchy(Parameters& rValues) {
    KRATOS_ERROR << "Calling base "
                    "FluidConstitutiveLaw::CalculateMaterialResponseCauchy "
                    "method. This class should not be instantiated. Please "
                    "check your constitutive law."
                 << std::endl;
}

int FluidConstitutiveLaw::Check(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_ERROR << "Calling base "
                    "FluidConstitutiveLaw::Check "
                    "method. This class should not be instantiated. Please "
                    "check your constitutive law."
                 << std::endl;
    return 999;
}

// Access /////////////////////////////////////////////////////////////////////

int& FluidConstitutiveLaw::CalculateValue(ConstitutiveLaw::Parameters& rParameters, const Variable<int>& rThisVariable,int& rValue) {
    return ConstitutiveLaw::CalculateValue(rParameters,rThisVariable,rValue);
}

double& FluidConstitutiveLaw::CalculateValue(ConstitutiveLaw::Parameters& rParameters, const Variable<double>& rThisVariable, double& rValue) {
    rValue = this->GetEffectiveViscosity(rParameters);
    return rValue;
}

Vector& FluidConstitutiveLaw::CalculateValue(ConstitutiveLaw::Parameters& rParameters, const Variable<Vector>& rThisVariable, Vector& rValue) {
    return ConstitutiveLaw::CalculateValue(rParameters,rThisVariable,rValue);
}

Matrix& FluidConstitutiveLaw::CalculateValue(ConstitutiveLaw::Parameters& rParameters, const Variable<Matrix>& rThisVariable, Matrix& rValue) {
    return ConstitutiveLaw::CalculateValue(rParameters,rThisVariable,rValue);
}

array_1d<double, 3 > & FluidConstitutiveLaw::CalculateValue(ConstitutiveLaw::Parameters& rParameters, const Variable<array_1d<double, 3 > >& rThisVariable,array_1d<double, 3 > & rValue) {
    return ConstitutiveLaw::CalculateValue(rParameters,rThisVariable,rValue);
}

array_1d<double, 6 > & FluidConstitutiveLaw::CalculateValue(ConstitutiveLaw::Parameters& rParameters, const Variable<array_1d<double, 6 > >& rThisVariable, array_1d<double, 6 > & rValue) {
    return ConstitutiveLaw::CalculateValue(rParameters,rThisVariable,rValue);
}

// Inquiry ////////////////////////////////////////////////////////////////////

ConstitutiveLaw::SizeType FluidConstitutiveLaw::WorkingSpaceDimension() {
    KRATOS_ERROR << "Calling base "
                    "FluidConstitutiveLaw::WorkingSpaceDimension "
                    "method. This class should not be instantiated. Please "
                    "check your constitutive law."
                 << std::endl;
    return 0;
}

ConstitutiveLaw::SizeType FluidConstitutiveLaw::GetStrainSize() const {
    KRATOS_ERROR << "Calling base "
                    "FluidConstitutiveLaw::GetStrainSize "
                    "method. This class should not be instantiated. Please "
                    "check your constitutive law."
                 << std::endl;
    return 0;
}

// Info ///////////////////////////////////////////////////////////////////////


std::string FluidConstitutiveLaw::Info() const {
    return "FluidConstitutiveLaw";
}

void FluidConstitutiveLaw::PrintInfo(std::ostream& rOStream) const {
    rOStream << this->Info();
}

void FluidConstitutiveLaw::PrintData(std::ostream& rOStream) const {
    rOStream << this->Info();
}

// Protected operations ///////////////////////////////////////////////////////

void FluidConstitutiveLaw::NewtonianConstitutiveMatrix2D(
    double EffectiveViscosity,
    Matrix& rC) {

    constexpr double two_thirds = 2./3.;
    constexpr double four_thirds = 4./3.;

    rC(0,0) = EffectiveViscosity * four_thirds;
    rC(0,1) = -EffectiveViscosity * two_thirds;
    rC(0,2) = 0.0;
    rC(1,0) = -EffectiveViscosity * two_thirds;
    rC(1,1) = EffectiveViscosity * four_thirds;
    rC(1,2) = 0.0;
    rC(2,0) = 0.0;
    rC(2,1) = 0.0;
    rC(2,2) = EffectiveViscosity;
}

void FluidConstitutiveLaw::NewtonianConstitutiveMatrix3D(
    double EffectiveViscosity,
    Matrix& rC) {

    rC.clear();

    constexpr double two_thirds = 2./3.;
    constexpr double four_thirds = 4./3.;

    rC(0,0) = EffectiveViscosity * four_thirds;
    rC(0,1) = -EffectiveViscosity * two_thirds;
    rC(0,2) = -EffectiveViscosity * two_thirds;

    rC(1,0) = -EffectiveViscosity * two_thirds;
    rC(1,1) = EffectiveViscosity * four_thirds;
    rC(1,2) = -EffectiveViscosity * two_thirds;

    rC(2,0) = -EffectiveViscosity * two_thirds;
    rC(2,1) = -EffectiveViscosity * two_thirds;
    rC(2,2) = EffectiveViscosity * four_thirds;

    rC(3,3) = EffectiveViscosity;
    rC(4,4) = EffectiveViscosity;
    rC(5,5) = EffectiveViscosity;
}

// Protected access ///////////////////////////////////////////////////////////

double FluidConstitutiveLaw::GetEffectiveViscosity(ConstitutiveLaw::Parameters& rParameters) const {
    KRATOS_ERROR << "Accessing base class FluidConstitutiveLaw::GetEffectiveViscosity." << std::endl;
    return 0.0;
}

double FluidConstitutiveLaw::GetValueFromTable(
    const Variable<double> &rVariableInput,
    const Variable<double> &rVariableOutput,
    ConstitutiveLaw::Parameters &rParameters) const
{
    // Get material properties from constitutive law parameters
    const Properties &r_properties = rParameters.GetMaterialProperties();

    double gauss_output;
    if (r_properties.HasTable(rVariableInput,rVariableOutput)) {
        // Get geometry and Gauss pt. data
        const auto &r_geom = rParameters.GetElementGeometry();
        const auto &r_N = rParameters.GetShapeFunctionsValues();

        // Compute the input variable Gauss pt. value
        double gauss_input = 0.0;
        for (unsigned int i_node = 0; i_node < r_N.size(); ++i_node) {
            const double &r_val = r_geom[i_node].FastGetSolutionStepValue(rVariableInput);
            gauss_input += r_val * r_N[i_node];
        }

        // Retrieve the output variable from the table
        const auto &r_table = r_properties.GetTable(rVariableInput, rVariableOutput);
        gauss_output = r_table.GetValue(gauss_input);
    } else {
        KRATOS_ERROR << "FluidConstitutiveLaw " << this->Info() << " has no table with variables " << rVariableInput << " " << rVariableOutput << std::endl;
    }

    return gauss_output;
}

// Serialization //////////////////////////////////////////////////////////////

void FluidConstitutiveLaw::save(Serializer& rSerializer) const {
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, ConstitutiveLaw );
}

void FluidConstitutiveLaw::load(Serializer& rSerializer) {
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, ConstitutiveLaw );
}

}