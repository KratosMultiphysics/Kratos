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
#include "thermal_elastic_isotropic_3d.h"
#include "custom_utilities/constitutive_law_utilities.h"
#include "custom_utilities/advanced_constitutive_law_utilities.h"

#include "includes/checks.h"
#include "structural_mechanics_application_variables.h"

namespace Kratos
{

/***********************************************************************************/
/***********************************************************************************/

void ThermalElasticIsotropic3D::CalculateMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{
    KRATOS_TRY

    Flags& r_constitutive_law_options = rValues.GetOptions();
    ConstitutiveLaw::StrainVectorType& r_strain_vector = rValues.GetStrainVector();

    if (r_constitutive_law_options.IsNot(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN)) {
        // Since we are in small strains, any strain measure works, e.g. CAUCHY_GREEN
        CalculateCauchyGreenStrain(rValues, r_strain_vector);
    }

    // We add the thermal contribution
    SubstractThermalStrain(r_strain_vector, mReferenceTemperature, rValues);

    // We add the initial strains
    AddInitialStrainVectorContribution<StrainVectorType>(r_strain_vector);


    if (r_constitutive_law_options.Is(ConstitutiveLaw::COMPUTE_STRESS)) {
        ConstitutiveLaw::StressVectorType &r_stress_vector = rValues.GetStressVector();
        CalculatePK2Stress(r_strain_vector, r_stress_vector, rValues);
        AddInitialStressVectorContribution<StressVectorType>(r_stress_vector);
    }

    if (r_constitutive_law_options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
        ConstitutiveLaw::VoigtSizeMatrixType &r_constitutive_matrix = rValues.GetConstitutiveMatrix();
        CalculateElasticMatrix(r_constitutive_matrix, rValues);
    }

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void ThermalElasticIsotropic3D::SubstractThermalStrain(
    ConstitutiveLaw::StrainVectorType &rStrainVector,
    const double ReferenceTemperature,
    ConstitutiveLaw::Parameters &rParameters,
    const bool IsPlaneStrain
    )
{
    AdvancedConstitutiveLawUtilities<6>::SubstractThermalStrain(rStrainVector, ReferenceTemperature, rParameters);
}

/***********************************************************************************/
/***********************************************************************************/

int ThermalElasticIsotropic3D::Check(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const ProcessInfo& rCurrentProcessInfo
    ) const
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

double& ThermalElasticIsotropic3D::CalculateValue(
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

void ThermalElasticIsotropic3D::InitializeMaterial(
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

/***********************************************************************************/
/***********************************************************************************/

void ThermalElasticIsotropic3D::CalculatePK2Stress(
    const Vector& rStrainVector,
    ConstitutiveLaw::StressVectorType& rStressVector,
    ConstitutiveLaw::Parameters& rValues
    )
{
    const Properties& r_material_properties = rValues.GetMaterialProperties();
    const auto &r_geom = rValues.GetElementGeometry();
    const auto &r_N = rValues.GetShapeFunctionsValues();
    const auto &r_process_info = rValues.GetProcessInfo();
    const double E  = r_material_properties.GetValue(YOUNG_MODULUS, r_geom, r_N, r_process_info);
    const double NU = r_material_properties.GetValue(POISSON_RATIO, r_geom, r_N, r_process_info);
    ConstitutiveLawUtilities<6>::CalculatePK2StressFromStrain(rStressVector, rStrainVector, E, NU);
}

/***********************************************************************************/
/***********************************************************************************/

void ThermalElasticIsotropic3D::CalculateElasticMatrix(
    ConstitutiveLaw::VoigtSizeMatrixType& rConstitutiveMatrix,
    ConstitutiveLaw::Parameters& rValues
    )
{
    const Properties& r_material_properties = rValues.GetMaterialProperties();
    const auto &r_geom = rValues.GetElementGeometry();
    const auto &r_N = rValues.GetShapeFunctionsValues();
    const auto &r_process_info = rValues.GetProcessInfo();
    const double E  = r_material_properties.GetValue(YOUNG_MODULUS, r_geom, r_N, r_process_info);
    const double NU = r_material_properties.GetValue(POISSON_RATIO, r_geom, r_N, r_process_info);
    ConstitutiveLawUtilities<6>::CalculateElasticMatrix(rConstitutiveMatrix, E, NU);
}

/***********************************************************************************/
/***********************************************************************************/

} // Namespace Kratos
