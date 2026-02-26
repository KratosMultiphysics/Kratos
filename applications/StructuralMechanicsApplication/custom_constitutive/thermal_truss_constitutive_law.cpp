// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    
//
// System includes

// External includes

// Project includes
#include "includes/checks.h"
#include "includes/properties.h"
#include "custom_constitutive/thermal_truss_constitutive_law.h"
#include "structural_mechanics_application_variables.h"

namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

ThermalTrussConstitutiveLaw::ThermalTrussConstitutiveLaw()
    : ConstitutiveLaw()
{
}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

ThermalTrussConstitutiveLaw::ThermalTrussConstitutiveLaw(const ThermalTrussConstitutiveLaw& rOther)
    : ConstitutiveLaw(rOther)
{
}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer ThermalTrussConstitutiveLaw::Clone() const
{
    return Kratos::make_shared<ThermalTrussConstitutiveLaw>(*this);
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

ThermalTrussConstitutiveLaw::~ThermalTrussConstitutiveLaw()
{
    // TODO: Add if necessary
}

//*************************CONSTITUTIVE LAW GENERAL FEATURES *************************
//************************************************************************************

void ThermalTrussConstitutiveLaw::GetLawFeatures(Features& rFeatures)
{
    //Set the strain size
    rFeatures.mStrainSize = 1;

    //Set the spacedimension
    rFeatures.mSpaceDimension = 3;
}
//************************************************************************************
//************************************************************************************

array_1d<double, 3 > & ThermalTrussConstitutiveLaw::GetValue(
    const Variable<array_1d<double, 3 > >& rThisVariable,
    array_1d<double, 3 > & rValue)
{
    KRATOS_ERROR << "Can't get the specified value" << std::endl;
    return rValue;
}

//************************************************************************************
//************************************************************************************

double& ThermalTrussConstitutiveLaw::CalculateValue(
    ConstitutiveLaw::Parameters& rParameterValues,
    const Variable<double>& rThisVariable,
    double& rValue)
{
    std::cout << "ThermalTrussConstitutiveLaw::CalculateValue double called" << std::endl;
    if(rThisVariable == TANGENT_MODULUS) 
    {
        std::cout << "ThermalTrussConstitutiveLaw::CalculateValue TANGENT_MODULUS called" << std::endl;
        rValue = rParameterValues.GetMaterialProperties()[YOUNG_MODULUS];
    }
    else if (rThisVariable == STRAIN_ENERGY){
        std::cout << "ThermalTrussConstitutiveLaw::CalculateValue STRAIN_ENERGY called" << std::endl;
        Vector current_strain = ZeroVector(1);
        rParameterValues.GetStrainVector(current_strain);
        rValue = 0.50 * rParameterValues.GetMaterialProperties()[YOUNG_MODULUS] * current_strain[0] * current_strain[0];
    }
    else KRATOS_ERROR << "Can't calculate the specified value" << std::endl;
    return rValue;
}

//************************************************************************************
//************************************************************************************

Vector& ThermalTrussConstitutiveLaw::CalculateValue(
    ConstitutiveLaw::Parameters& rParameterValues,
    const Variable<Vector>& rThisVariable,
    Vector& rValue)
{
    if(rThisVariable == NORMAL_STRESS)
    {
        std::cout << "ThermalTrussConstitutiveLaw::CalculateValue NORMAL_STRESS called" << std::endl;
        std::cout << "Hi normal stress called" << std::endl;
        const double current_stress = this->CalculateStressElastic(rParameterValues);
        constexpr SizeType dofs = 6;
        rValue = ZeroVector(dofs);
        rValue[0] = -1.0 * current_stress;
        rValue[3] = 1.0 * current_stress;
    }
    else KRATOS_ERROR << "Can't calculate the specified value" << std::endl;
    return rValue;
}

//************************************************************************************
//************************************************************************************

array_1d<double, 3 > & ThermalTrussConstitutiveLaw::CalculateValue(
    ConstitutiveLaw::Parameters& rParameterValues,
    const Variable<array_1d<double, 3 > >& rVariable,
	array_1d<double, 3 > & rValue)
    {
        if (rVariable == FORCE)
        {
            std::cout << "ThermalTrussConstitutiveLaw::CalculateValue FORCE called" << std::endl;
            std::cout << "Hi force called" << std::endl;
            constexpr SizeType dimension = 3;
            rValue = ZeroVector(dimension);
            //rValue[0] = this->mStressState;
            rValue[0] = this->CalculateStressElastic(rParameterValues);
            rValue[1] = 0.0;
            rValue[2] = 0.0;
        }
        else KRATOS_ERROR << "Can't calculate the specified value" << std::endl;
        return rValue;
    }

//************************************************************************************
//************************************************************************************
void ThermalTrussConstitutiveLaw::CalculateMaterialResponsePK2(Parameters& rValues)
{
    std::cout << "ThermalTrussConstitutiveLaw::CalculateMaterialResponsePK2 called" << std::endl;
    //std::cout << "Hi Mat response called" << std::endl;
    Vector& stress_vector = rValues.GetStressVector();
    if (stress_vector.size() != 1) stress_vector.resize(1, false);
    stress_vector[0] = this->CalculateStressElastic(rValues);
}
//************************************************************************************
//************************************************************************************

double ThermalTrussConstitutiveLaw::CalculateStressElastic(
    ConstitutiveLaw::Parameters& rParameterValues)
{
    std::cout << "ThermalTrussConstitutiveLaw::CalculateStressElastic called" << std::endl;
    //std::cout << "Hi stress elastic called" << std::endl;
    Vector current_strain = ZeroVector(1);
    rParameterValues.GetStrainVector(current_strain);

    // Add thermal contribution
    double alpha = rParameterValues.GetMaterialProperties()[THERMAL_EXPANSION_COEFFICIENT];
    double T_ref = rParameterValues.GetMaterialProperties()[REFERENCE_TEMPERATURE];
    const auto& r_geometry = rParameterValues.GetElementGeometry();
    double T_node_1 = r_geometry[0].GetSolutionStepValue(TEMPERATURE);
    double T_node_2 = r_geometry[1].GetSolutionStepValue(TEMPERATURE);
    const double current_temperature_gp = (T_node_1 + T_node_2)* 0.5;
    alpha *= (current_temperature_gp - T_ref);
    //std::cout << "thermal strian is " << alpha << std::endl;
    current_strain[0] -= 1.0*alpha;

    double tangent_modulus(0.0);
    CalculateValue(rParameterValues,TANGENT_MODULUS,tangent_modulus);

    const double current_stress = tangent_modulus*current_strain[0];
    //std::cout << " stress is " << current_stress << std::endl;
    return current_stress;
}

//************************************************************************************
//************************************************************************************

int ThermalTrussConstitutiveLaw::Check(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const ProcessInfo& rCurrentProcessInfo
) const
{
    KRATOS_ERROR_IF_NOT(rElementGeometry[0].SolutionStepsDataHas(TEMPERATURE))  << "The TEMPERATURE variable is not available at the nodes." << std::endl;
    KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(THERMAL_EXPANSION_COEFFICIENT)) << "The THERMAL_EXPANSION_COEFFICIENT is not set in the material properties." << std::endl;
    KRATOS_ERROR_IF(rMaterialProperties[THERMAL_EXPANSION_COEFFICIENT] < 0.0)   << "The THERMAL_EXPANSION_COEFFICIENT is negative..." << std::endl;
    KRATOS_ERROR_IF_NOT(rElementGeometry.Has(REFERENCE_TEMPERATURE) || rMaterialProperties.Has(REFERENCE_TEMPERATURE)) << "The REFERENCE_TEMPERATURE is not given in the material properties nor via SetValue()" << std::endl;
    KRATOS_CHECK(rMaterialProperties.Has(YOUNG_MODULUS));
    KRATOS_ERROR_IF(YOUNG_MODULUS.Key() == 0 || rMaterialProperties[YOUNG_MODULUS] <= 0.00)
     << "YOUNG_MODULUS has Key zero or invalid value " << std::endl;

    KRATOS_ERROR_IF(DENSITY.Key() == 0 || rMaterialProperties[DENSITY] < 0.00)
     << "DENSITY has Key zero or invalid value " << std::endl;

    return 0;

}

} // Namespace Kratos
