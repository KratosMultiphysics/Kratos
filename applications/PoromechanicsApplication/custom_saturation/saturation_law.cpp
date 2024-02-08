//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ignasi de Pouplana
//                   Danilo Cavalcanti
//

/* Project includes */
#include "custom_saturation/saturation_law.hpp"

namespace Kratos
{

SaturationLaw::Pointer SaturationLaw::Clone() const
{
    KRATOS_ERROR <<  "Called the virtual function for Clone"<< std::endl;
}

//------------------------------------------------------------------------------------------------

// SaturationLaw::SizeType SaturationLaw::WorkingSpaceDimension()
// {
//     KRATOS_ERROR <<  "Called the virtual function for WorkingSpaceDimension"<< std::endl;
// }

//------------------------------------------------------------------------------------------------

bool SaturationLaw::Has(const Variable<bool>& rThisVariable)
{
    return false;
}
bool SaturationLaw::Has(const Variable<int>& rThisVariable)
{
    return false;
}
bool SaturationLaw::Has(const Variable<double>& rThisVariable)
{
    return false;
}
bool SaturationLaw::Has(const Variable<Vector>& rThisVariable)
{
    return false;
}
bool SaturationLaw::Has(const Variable<Matrix>& rThisVariable)
{
    return false;
}
bool SaturationLaw::Has(const Variable<array_1d<double, 3 > >& rThisVariable)
{
    return false;
}
bool SaturationLaw::Has(const Variable<array_1d<double, 6 > >& rThisVariable)
{
    return false;
}

//------------------------------------------------------------------------------------------------

bool& SaturationLaw::GetValue(const Variable<bool>& rThisVariable, bool& rValue)
{
    return rValue;
}
int& SaturationLaw::GetValue(const Variable<int>& rThisVariable, int& rValue)
{
    return rValue;
}
double& SaturationLaw::GetValue(const Variable<double>& rThisVariable, double& rValue)
{
    return rValue;
}
Vector& SaturationLaw::GetValue(const Variable<Vector>& rThisVariable, Vector& rValue)
{
    return rValue;
}
Matrix& SaturationLaw::GetValue(const Variable<Matrix>& rThisVariable, Matrix& rValue)
{
    return rValue;
}
array_1d<double, 3 > & SaturationLaw::GetValue(const Variable<array_1d<double, 3 > >& rThisVariable,
        array_1d<double, 3 > & rValue)
{
    return rValue;
}
array_1d<double, 6 > & SaturationLaw::GetValue(const Variable<array_1d<double, 6 > >& rThisVariable,
        array_1d<double, 6 > & rValue)
{
    return rValue;
}

//------------------------------------------------------------------------------------------------

void SaturationLaw::SetValue(const Variable<bool>& rThisVariable,
                               const bool& Value,
                               const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_ERROR <<  "Called the virtual function for SetValue"<< std::endl;
}
void SaturationLaw::SetValue(const Variable<int>& rThisVariable,
                               const int& Value,
                               const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_ERROR <<  "Called the virtual function for SetValue"<< std::endl;
}
void SaturationLaw::SetValue(const Variable<double>& rVariable,
                               const double& rValue,
                               const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_ERROR <<  "Called the virtual function for SetValue"<< std::endl;
}
void SaturationLaw::SetValue(const Variable<Vector >& rVariable,
                               const Vector& rValue, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_ERROR <<  "Called the virtual function for SetValue"<< std::endl;
}
void SaturationLaw::SetValue(const Variable<Matrix >& rVariable,
                               const Matrix& rValue, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_ERROR <<  "Called the virtual function for SetValue"<< std::endl;
}
void SaturationLaw::SetValue(const Variable<array_1d<double, 3 > >& rVariable,
                               const array_1d<double, 3 > & rValue,
                               const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_ERROR <<  "Called the virtual function for SetValue"<< std::endl;
}
void SaturationLaw::SetValue(const Variable<array_1d<double, 6 > >& rVariable,
                               const array_1d<double, 6 > & rValue,
                               const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_ERROR <<  "Called the virtual function for SetValue"<< std::endl;
}

//------------------------------------------------------------------------------------------------

bool& SaturationLaw::CalculateValue(Parameters& rParameterValues, const Variable<bool>& rThisVariable, bool& rValue)
{
    return rValue;
}
int& SaturationLaw::CalculateValue(Parameters& rParameterValues, const Variable<int>& rThisVariable, int& rValue)
{
    return rValue;
}
double& SaturationLaw::CalculateValue(Parameters& rParameterValues, const Variable<double>& rThisVariable, double& rValue)
{
    return rValue;
}
Vector& SaturationLaw::CalculateValue(Parameters& rParameterValues, const Variable<Vector>& rThisVariable, Vector& rValue)
{
    return rValue;
}
Matrix& SaturationLaw::CalculateValue(Parameters& rParameterValues, const Variable<Matrix>& rThisVariable, Matrix& rValue)
{
    return rValue;
}
array_1d<double, 3 > & SaturationLaw::CalculateValue(Parameters& rParameterValues, const Variable<array_1d<double, 3 > >& rVariable,
        array_1d<double, 3 > & rValue)
{
    return rValue;
}
array_1d<double, 6 > & SaturationLaw::CalculateValue(Parameters& rParameterValues, const Variable<array_1d<double, 6 > >& rVariable,
        array_1d<double, 6 > & rValue)
{
    return rValue;
}

//------------------------------------------------------------------------------------------------

int SaturationLaw::Check(const Properties& rMaterialProperties,
                           const GeometryType& rElementGeometry,
                           const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    // Verify Properties variables
    if(rMaterialProperties.Has(RESIDUAL_WATER_SATURATION)) {
        const double& residual_water_saturation = rMaterialProperties[RESIDUAL_WATER_SATURATION];
        const bool check = static_cast<bool>((residual_water_saturation < 0.0) || (residual_water_saturation > 1.0));
        KRATOS_ERROR_IF(check) << "RESIDUAL_WATER_SATURATION has an invalid value " << std::endl;
    } else {
        KRATOS_ERROR << "RESIDUAL_WATER_SATURATION not defined" << std::endl;
    }

    if(rMaterialProperties.Has(PORE_SIZE_FACTOR)) {
        KRATOS_ERROR_IF(rMaterialProperties[PORE_SIZE_FACTOR] <= 0.0) << "PORE_SIZE_FACTOR has an invalid value " << std::endl;
    } else {
        KRATOS_ERROR << "PORE_SIZE_FACTOR not defined" << std::endl;
    }

    if(!(rMaterialProperties.Has(GAS_ENTRY_PRESSURE))) {
        KRATOS_ERROR << "GAS_ENTRY_PRESSURE not defined" << std::endl;
    }

    return 0;

    KRATOS_CATCH("");
}

//------------------------------------------------------------------------------------------------

bool SaturationLaw::ValidateInput(const Properties& rMaterialProperties)
{
  return false;
}

//------------------------------------------------------------------------------------------------

void SaturationLaw::InitializeMaterial(const Properties& rMaterialProperties,
        const GeometryType& rElementGeometry,
        const Vector& rShapeFunctionsValues)
{
}

//------------------------------------------------------------------------------------------------

void SaturationLaw::CalculateMaterialResponse (Parameters& rValues)
{
    //Check
    rValues.CheckAllParameters();

    //Initialize main variables
    SaturationLawVariables Variables;
    this->InitializeSaturationLawVariables(Variables,rValues);

    this->CalculateWaterSaturationDegree(Variables,rValues);

    this->CalculateEffectiveSaturation(Variables,rValues);
    this->CalculateWaterRelativePermeability(Variables,rValues);
    this->CalculateGasRelativePermeability(Variables,rValues);
}

//------------------------------------------------------------------------------------------------

void SaturationLaw::InitializeMaterialResponse (Parameters& rValues)
{
    KRATOS_ERROR_IF(this->RequiresInitializeMaterialResponse()) <<  "Calling virtual function for InitializeMaterialResponse. Please implement InitializeMaterialResponse or RequiresInitializeMaterialResponse in case this CL does not require it" << std::endl;
}

//------------------------------------------------------------------------------------------------

void SaturationLaw::FinalizeMaterialResponse (Parameters& rValues)
{
    KRATOS_ERROR_IF(this->RequiresFinalizeMaterialResponse()) <<  "Calling virtual function for FinalizeMaterialResponse. Please implement FinalizeMaterialResponse or RequiresFinalizeMaterialResponse in case this CL does not require it" << std::endl;
}

//------------------------------------------------------------------------------------------------

void SaturationLaw::CalculateSaturation (Parameters& rValues)
{
    //Check
    rValues.CheckAllParameters();

    //Initialize main variables
    SaturationLawVariables Variables;
    this->InitializeSaturationLawVariables(Variables,rValues);

    this->CalculateWaterSaturationDegree(Variables,rValues);
}

//------------------------------------------------------------------------------------------------

void SaturationLaw::ResetMaterial(const Properties& rMaterialProperties,
                                    const GeometryType& rElementGeometry,
                                    const Vector& rShapeFunctionsValues)
{
    KRATOS_ERROR <<  "Calling virtual function for ResetMaterial"<< std::endl;
}

//------------------------------------------------------------------------------------------------

void SaturationLaw::InitializeSaturationLawVariables(SaturationLawVariables& rVariables, Parameters& rValues)
{
    const Properties& properties = rValues.GetMaterialProperties();
    const Vector& N = rValues.GetShapeFunctionsValues();
    const Element::GeometryType& geometry = rValues.GetElementGeometry();
    const unsigned int number_of_nodes = geometry.size();

    rVariables.Swr = properties[RESIDUAL_WATER_SATURATION];
    rVariables.lambda = properties[PORE_SIZE_FACTOR];
    rVariables.pb = properties[GAS_ENTRY_PRESSURE];
    // Capilar pore pressure: pc = pg - pw
    rVariables.pc = 0.0;
    for (unsigned int i = 0; i < number_of_nodes; i++) {
        rVariables.pc += N[i] * (geometry[i].GetSolutionStepValue(GAS_PRESSURE) - geometry[i].GetSolutionStepValue(WATER_PRESSURE));
    }
}

//------------------------------------------------------------------------------------------------

void SaturationLaw::CalculateWaterSaturationDegree (SaturationLawVariables& rVariables, Parameters& rValues)
{
    KRATOS_ERROR <<  "Calling virtual function for CalculateWaterSaturationDegree"<< std::endl;
}

//------------------------------------------------------------------------------------------------

void SaturationLaw::CalculateEffectiveSaturation(SaturationLawVariables& rVariables, Parameters& rValues)
{
    const double& Sw = rValues.GetSw();
    rVariables.Se = (Sw - rVariables.Swr)/(1.0 - rVariables.Swr);
}

//------------------------------------------------------------------------------------------------

void SaturationLaw::CalculateWaterRelativePermeability (SaturationLawVariables& rVariables, Parameters& rValues)
{
    KRATOS_ERROR <<  "Calling virtual function for CalculateWaterRelativePermeability"<< std::endl;
}

//------------------------------------------------------------------------------------------------

void SaturationLaw::CalculateGasRelativePermeability (SaturationLawVariables& rVariables, Parameters& rValues)
{
    KRATOS_ERROR <<  "Calling virtual function for CalculateGasRelativePermeability"<< std::endl;
}

//------------------------------------------------------------------------------------------------

void SaturationLaw::save(Serializer& rSerializer) const
{
    // there are no member variables to be saved
}

void SaturationLaw::load(Serializer& rSerializer)
{
    // there are no member variables to be loaded
}

}