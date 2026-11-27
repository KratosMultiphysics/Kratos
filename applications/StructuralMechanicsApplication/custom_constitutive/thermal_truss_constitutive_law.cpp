// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Talhah Ansari
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
    : TrussConstitutiveLaw()
{
}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

ThermalTrussConstitutiveLaw::ThermalTrussConstitutiveLaw(const ThermalTrussConstitutiveLaw& rOther)
    : TrussConstitutiveLaw(rOther),
      mReferenceTemperature(rOther.mReferenceTemperature)
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

double ThermalTrussConstitutiveLaw::CalculateStressElastic(
    ConstitutiveLaw::Parameters& rParameterValues)
{
    Vector current_strain(1);
    rParameterValues.GetStrainVector(current_strain);

    // Add thermal contribution
    const auto& r_props = rParameterValues.GetMaterialProperties();
    const double alpha = r_props[THERMAL_EXPANSION_COEFFICIENT];
    const auto& r_geometry = rParameterValues.GetElementGeometry();
    const double current_temperature_gp = 0.5 * (
        r_geometry[0].FastGetSolutionStepValue(TEMPERATURE) +
        r_geometry[1].FastGetSolutionStepValue(TEMPERATURE));
    const double thermal_strain = alpha * (current_temperature_gp - mReferenceTemperature);
    current_strain[0] -= 1.0*thermal_strain;

    double tangent_modulus(0.0);
    CalculateValue(rParameterValues,TANGENT_MODULUS,tangent_modulus);

    const double current_stress = tangent_modulus*current_strain[0];

    if (rParameterValues.GetOptions().Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
        Matrix& r_const_matrix = rParameterValues.GetConstitutiveMatrix();
        r_const_matrix.resize(1, 1, false);
        r_const_matrix(0, 0) = tangent_modulus;
    }
    return current_stress;
}

void ThermalTrussConstitutiveLaw::InitializeMaterial(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const Vector&)
{
    if (rElementGeometry.Has(REFERENCE_TEMPERATURE)) {
        mReferenceTemperature = rElementGeometry.GetValue(REFERENCE_TEMPERATURE);
    } else if (rMaterialProperties.Has(REFERENCE_TEMPERATURE)) {
        mReferenceTemperature = rMaterialProperties[REFERENCE_TEMPERATURE];
    }
}

//************************************************************************************
//************************************************************************************

int ThermalTrussConstitutiveLaw::Check(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const ProcessInfo& rCurrentProcessInfo
) const
{
    TrussConstitutiveLaw::Check(rMaterialProperties, rElementGeometry, rCurrentProcessInfo);

    KRATOS_ERROR_IF(rElementGeometry.size() != 2)
        << "ThermalTrussConstitutiveLaw only supports two-node geometries." << std::endl;
    KRATOS_ERROR_IF(rMaterialProperties.Has(INTEGRATION_ORDER) && rMaterialProperties[INTEGRATION_ORDER] != 1)
        << "ThermalTrussConstitutiveLaw requires one integration point; set INTEGRATION_ORDER to 1." << std::endl;

    KRATOS_ERROR_IF_NOT(rElementGeometry[0].SolutionStepsDataHas(TEMPERATURE))  << "The TEMPERATURE variable is not available at the nodes." << std::endl;
    KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(THERMAL_EXPANSION_COEFFICIENT)) << "The THERMAL_EXPANSION_COEFFICIENT is not set in the material properties." << std::endl;
    KRATOS_ERROR_IF(rMaterialProperties[THERMAL_EXPANSION_COEFFICIENT] < 0.0)   << "The THERMAL_EXPANSION_COEFFICIENT is negative..." << std::endl;
    KRATOS_ERROR_IF_NOT(rElementGeometry.Has(REFERENCE_TEMPERATURE) || rMaterialProperties.Has(REFERENCE_TEMPERATURE)) << "The REFERENCE_TEMPERATURE is not given in the material properties nor via SetValue()" << std::endl;

    return 0;

}

} // Namespace Kratos
