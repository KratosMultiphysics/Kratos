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
#include "custom_constitutive_thermal/thermal_elastic_isotropic_3d.h"

namespace Kratos
{

/***********************************************************************************/
/***********************************************************************************/

void  ThermalElasticIsotropic3D::CalculateMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{
    KRATOS_TRY;
    const auto& r_process_info = rValues.GetProcessInfo();
    if (r_process_info[STEP] == 1) // We storage the ref temperature as the initial one
        mReferenceTemperature = AdvancedConstitutiveLawUtilities<6>::CalculateInGaussPoint(TEMPERATURE, rValues);

    Flags& r_constitutive_law_options = rValues.GetOptions();
    ConstitutiveLaw::StrainVectorType& r_strain_vector = rValues.GetStrainVector();

    const bool use_elem_provided_flag_backup = r_constitutive_law_options.Is(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);

    if (r_constitutive_law_options.IsNot(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN)) {
        // Since we are in small strains, any strain measure works, e.g. CAUCHY_GREEN
        CalculateCauchyGreenStrain(rValues, r_strain_vector);
    }

    // We add the thermal contribution
    AdvancedConstitutiveLawUtilities<6>::SubstractThermalStrain(r_strain_vector, mReferenceTemperature, rValues);

    // We force to use the already computed strain in the base CL
    r_constitutive_law_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED, true);

    BaseType::CalculateMaterialResponsePK2(rValues);

    // We reset the flag
    r_constitutive_law_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED, use_elem_provided_flag_backup);

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

int ThermalElasticIsotropic3D::Check(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_ERROR_IF_NOT(rElementGeometry[0].HasNodalSolutionStepVariable(TEMPERATURE)) << "The TEMPERATURE variable is not available at the nodes." << std::endl;
    KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(THERMAL_EXPANSION_COEFFICIENT))        << "The THERMAL_EXPANSION_COEFFICIENT is not set in the material properties." << std::endl;
    KRATOS_ERROR_IF(rMaterialProperties[THERMAL_EXPANSION_COEFFICIENT] < 0.0)          << "The THERMAL_EXPANSION_COEFFICIENT is negative..." << std::endl;
    BaseType::Check(rMaterialProperties, rElementGeometry, rCurrentProcessInfo);
    return 0;
}

/***********************************************************************************/
/***********************************************************************************/

} // Namespace Kratos
