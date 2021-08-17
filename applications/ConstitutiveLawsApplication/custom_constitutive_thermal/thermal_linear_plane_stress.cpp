// KRATOS ___                _   _ _         _   _             __                       _
//       / __\___  _ __  ___| |_(_) |_ _   _| |_(_)_   _____  / /  __ ___      _____   /_\  _ __  _ __
//      / /  / _ \| '_ \/ __| __| | __| | | | __| \ \ / / _ \/ /  / _` \ \ /\ / / __| //_\\| '_ \| '_  |
//     / /__| (_) | | | \__ \ |_| | |_| |_| | |_| |\ V /  __/ /__| (_| |\ V  V /\__ \/  _  \ |_) | |_) |
//     \____/\___/|_| |_|___/\__|_|\__|\__,_|\__|_| \_/ \___\____/\__,_| \_/\_/ |___/\_/ \_/ .__/| .__/
//                                                                                         |_|   |_|
//
//  License:         BSD License
//                   license: structural_mechanics_application/license.txt
//
//  Main authors:    Alejandro Cornejo
//

// System includes
#include <iostream>

// External includes

// Project includes
#include "custom_constitutive_thermal/thermal_linear_plane_stress.h"

#include "includes/checks.h"
#include "structural_mechanics_application_variables.h"

namespace Kratos
{

/***********************************************************************************/
/***********************************************************************************/

void ThermalLinearPlaneStress::CalculateMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{
    KRATOS_TRY

    Flags& r_constitutive_law_options = rValues.GetOptions();
    ConstitutiveLaw::StrainVectorType& r_strain_vector = rValues.GetStrainVector();

    const bool use_elem_provided_flag_backup = r_constitutive_law_options.Is(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);

    if (r_constitutive_law_options.IsNot(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN)) {
        // Since we are in small strains, any strain measure works, e.g. CAUCHY_GREEN
        CalculateCauchyGreenStrain(rValues, r_strain_vector);
    }

    // We add the thermal contribution
    AdvancedConstitutiveLawUtilities<3>::SubstractThermalStrain(r_strain_vector, mReferenceTemperature, rValues);

    // We force to use the already computed strain in the base CL
    r_constitutive_law_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);

    BaseType::CalculateMaterialResponsePK2(rValues);

    // We reset the flag
    r_constitutive_law_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, use_elem_provided_flag_backup);

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

int ThermalLinearPlaneStress::Check(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_ERROR_IF_NOT(rElementGeometry[0].SolutionStepsDataHas(TEMPERATURE))  << "The TEMPERATURE variable is not available at the nodes." << std::endl;
    KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(THERMAL_EXPANSION_COEFFICIENT)) << "The THERMAL_EXPANSION_COEFFICIENT is not set in the material properties." << std::endl;
    KRATOS_ERROR_IF(rMaterialProperties[THERMAL_EXPANSION_COEFFICIENT] < 0.0)   << "The THERMAL_EXPANSION_COEFFICIENT is negative..." << std::endl;
    KRATOS_ERROR_IF_NOT(rElementGeometry.Has(REFERENCE_TEMPERATURE) || rMaterialProperties.Has(REFERENCE_TEMPERATURE)) << "The REFERENCE_TEMPERATURE is not given in the material properties nor via SetValue()" << std::endl;
    BaseType::Check(rMaterialProperties, rElementGeometry, rCurrentProcessInfo);
    return 0;
}

/***********************************************************************************/
/***********************************************************************************/

Vector& ThermalLinearPlaneStress::CalculateValue(
    ConstitutiveLaw::Parameters& rParameterValues,
    const Variable<Vector>& rThisVariable,
    Vector& rValue
    )
{
    if (rThisVariable == STRAIN ||
        rThisVariable == GREEN_LAGRANGE_STRAIN_VECTOR ||
        rThisVariable == ALMANSI_STRAIN_VECTOR) {

        Flags& r_flags = rParameterValues.GetOptions();

        // Previous flags saved
        const bool flag_const_tensor = r_flags.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
        const bool flag_stress = r_flags.Is(ConstitutiveLaw::COMPUTE_STRESS);

        // Set flags to only compute the strain
        r_flags.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);
        r_flags.Set(ConstitutiveLaw::COMPUTE_STRESS, false);

        ThermalLinearPlaneStress::CalculateMaterialResponsePK2(rParameterValues);
        if (rValue.size() != GetStrainSize()) {
            rValue.resize(GetStrainSize());
        }
        noalias(rValue) = rParameterValues.GetStrainVector();

        // Previous flags restored
        r_flags.Set( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, flag_const_tensor );
        r_flags.Set( ConstitutiveLaw::COMPUTE_STRESS, flag_stress );

    } else if (rThisVariable == STRESSES ||
        rThisVariable == CAUCHY_STRESS_VECTOR ||
        rThisVariable == KIRCHHOFF_STRESS_VECTOR ||
        rThisVariable == PK2_STRESS_VECTOR) {

        Flags& r_flags = rParameterValues.GetOptions();

        // Previous flags saved
        const bool flag_const_tensor = r_flags.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
        const bool flag_stress = r_flags.Is(ConstitutiveLaw::COMPUTE_STRESS);

        // Set flags to only compute the stress
        r_flags.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);
        r_flags.Set(ConstitutiveLaw::COMPUTE_STRESS, true);

        ThermalLinearPlaneStress::CalculateMaterialResponsePK2(rParameterValues);
        if (rValue.size() != GetStrainSize()) {
            rValue.resize(GetStrainSize());
        }
        noalias(rValue) = rParameterValues.GetStressVector();

        // Previous flags restored
        r_flags.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, flag_const_tensor);
        r_flags.Set(ConstitutiveLaw::COMPUTE_STRESS, flag_stress);
    } else if (rThisVariable == INITIAL_STRAIN_VECTOR) {
        if (this->HasInitialState()) {
            if (rValue.size() != GetStrainSize()) {
                rValue.resize(GetStrainSize());
            }
            noalias(rValue) = GetInitialState().GetInitialStrainVector();
        } else {
            noalias(rValue) = ZeroVector(0);
        }
    }
    return (rValue);
}

/***********************************************************************************/
/***********************************************************************************/

double& ThermalLinearPlaneStress::CalculateValue(
    ConstitutiveLaw::Parameters& rParameterValues,
    const Variable<double>& rThisVariable, double& rValue
    )
{
    if (rThisVariable == REFERENCE_TEMPERATURE) {
        rValue = mReferenceTemperature;
    } else {
        BaseType::CalculateValue(rParameterValues, rThisVariable, rValue);
    }
    return (rValue);
}

/***********************************************************************************/
/***********************************************************************************/

void ThermalLinearPlaneStress::InitializeMaterial(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const Vector& rShapeFunctionsValues
    )
{
    if (rElementGeometry.Has(REFERENCE_TEMPERATURE)) {
        mReferenceTemperature = rElementGeometry.GetValue(REFERENCE_TEMPERATURE);
    } else if (rMaterialProperties.Has(REFERENCE_TEMPERATURE)) {
        mReferenceTemperature = rMaterialProperties[REFERENCE_TEMPERATURE];
    }
}

} // Namespace Kratos
